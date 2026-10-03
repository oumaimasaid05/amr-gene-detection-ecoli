import glob, os
import pandas as pd

cols = ["qseqid", "contig", "pident", "length", "qlen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore"]
names = pd.read_csv("data/strains.csv").set_index("accession")["strain"].to_dict()
amr_all = pd.read_csv("results/amr_hits_classified.csv")
amr_genes = amr_all[amr_all["Subtype"] == "AMR"]   # point mutations can't be found by BLAST of genes

def overlaps(c, lo, hi, df):
    d = df[df["Contig id"] == c]
    for _, r in d.iterrows():
        a, b = min(r["Start"], r["Stop"]), max(r["Start"], r["Stop"])
        if min(b, hi) - max(a, lo) > 0:
            return True
    return False

summary, missed, card_only = [], [], []
for f in sorted(glob.glob("results/blast/*.tsv")):
    acc = os.path.basename(f)[:-4]
    if os.path.getsize(f) == 0:
        continue
    b = pd.read_csv(f, sep="\t", names=cols)
    b["qcov"] = b["length"] / b["qlen"]
    b = b[(b["pident"] >= 90) & (b["qcov"] >= 0.8)].copy()
    b["lo"] = b[["sstart", "send"]].min(axis=1)
    b["hi"] = b[["sstart", "send"]].max(axis=1)
    b = b.sort_values("bitscore", ascending=False)
    kept = []
    for _, r in b.iterrows():
        dup = False
        for k in kept:
            if k["contig"] == r["contig"]:
                inter = min(k["hi"], r["hi"]) - max(k["lo"], r["lo"])
                if inter > 0.5 * min(k["hi"] - k["lo"], r["hi"] - r["lo"]):
                    dup = True
                    break
        if not dup:
            kept.append(r)
    card = pd.DataFrame(kept)
    card["card_gene"] = card["qseqid"].str.split("|").str[-1].str.split(" \\[").str[0]

    a_genes = amr_genes[amr_genes["accession"] == acc]
    a_all = amr_all[amr_all["accession"] == acc]
    found = 0
    for _, g in a_genes.iterrows():
        lo, hi = min(g["Start"], g["Stop"]), max(g["Start"], g["Stop"])
        if overlaps(g["Contig id"], lo, hi, card.rename(columns={"contig": "Contig id", "lo": "Start", "hi": "Stop"})):
            found += 1
        else:
            missed.append((names[acc], g["Element symbol"], g["category"]))
    for _, r in card.iterrows():
        if not overlaps(r["contig"], r["lo"], r["hi"], a_all):
            card_only.append((names[acc], r["card_gene"], round(r["pident"], 1), round(r["qcov"], 2)))
    summary.append({"strain": names[acc], "AMRFinder genes": len(a_genes),
                    "also found by CARD-BLAST": found,
                    "CARD-BLAST only": sum(1 for x in card_only if x[0] == names[acc])})

s = pd.DataFrame(summary)
s["agreement_%"] = (100 * s["also found by CARD-BLAST"] / s["AMRFinder genes"]).round(1)
s.to_csv("results/tool_concordance.csv", index=False)
print(s.to_string(index=False))
print("\nAMRFinderPlus genes NOT found by CARD-BLAST:")
for m in missed:
    print("  ", m)
co = pd.DataFrame(card_only, columns=["strain", "card_gene", "pident", "qcov"])
co.to_csv("results/card_only_hits.csv", index=False)
print("\nMost frequent CARD-BLAST-only genes:")
print(co["card_gene"].value_counts().head(15).to_string())
