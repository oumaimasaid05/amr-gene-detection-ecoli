import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

c = pd.read_csv("results/contigs.tsv", sep="\t", header=None,
                names=["accession", "header", "length"])
c["contig"] = c["header"].str.split().str[0]
longest = c.groupby("accession")["length"].transform("max")
c["location"] = (c["length"] == longest).map({True: "chromosome", False: "plasmid/other"})
c.to_csv("results/contigs_annotated.csv", index=False)

a = pd.read_csv("results/amr_hits_classified.csv")
a = a[(a["category"] == "acquired")]
a = a.merge(c[["accession", "contig", "length", "location"]],
            left_on=["accession", "Contig id"], right_on=["accession", "contig"], how="left")
a.to_csv("results/gene_locations.csv", index=False)
print("hits with no matching contig:", a["location"].isna().sum())

tab = (a.groupby(["strain", "location"]).size().unstack(fill_value=0)
         .reindex(columns=["chromosome", "plasmid/other"], fill_value=0))
tab = tab.loc[tab.sum(axis=1).sort_values().index]
print(tab)
tab.plot(kind="barh", stacked=True, figsize=(8, 5), color=["#718096", "#c53030"])
plt.xlabel("Acquired AMR genes")
plt.title("Acquired AMR genes: chromosome vs plasmid")
plt.xticks(range(0, int(tab.sum(axis=1).max()) + 2, 2))
plt.tight_layout()
plt.savefig("results/figures/chromosome_vs_plasmid.png", dpi=300)

print()
for (s, ctg), d in a.sort_values("Start").groupby(["strain", "Contig id"]):
    print(f"{s} | {ctg} | {d['length'].iloc[0]:,} bp | {d['location'].iloc[0]}")
    print("   " + "  ".join(f"{r['Element symbol']}@{int(r['Start'])}" for _, r in d.iterrows()))
