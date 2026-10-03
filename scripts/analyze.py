import glob, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

os.makedirs("results/figures", exist_ok=True)
names = pd.read_csv("data/strains.csv").set_index("accession")["strain"].to_dict()
strains = [names[a] for a in sorted(names)]

frames = []
for f in sorted(glob.glob("results/amrfinder/*.tsv")):
    d = pd.read_csv(f, sep="\t")
    d["accession"] = os.path.basename(f)[:-4]
    frames.append(d)
amr = pd.concat(frames, ignore_index=True)
amr = amr[amr["Type"] == "AMR"].copy()
amr["strain"] = amr["accession"].map(names)

INTRINSIC = {"acrF", "emrD", "mdtM"}
def category(r):
    s = r["Element symbol"]
    if s.startswith("blaEC") or s in INTRINSIC:
        return "intrinsic"
    if "POINT" in str(r["Subtype"]):
        return "mutation"
    return "acquired"
amr["category"] = amr.apply(category, axis=1)
amr.to_csv("results/amr_hits_classified.csv", index=False)

core = amr[amr["category"] != "intrinsic"].drop_duplicates(["strain", "Element symbol"]).copy()

# 1. presence/absence heatmap (acquired genes + mutations)
pres = (core.assign(v=1)
            .pivot_table(index="Element symbol", columns="strain", values="v",
                         aggfunc="max", fill_value=0)
            .reindex(columns=strains, fill_value=0))
pres.to_csv("results/presence_absence_matrix.csv")
g = sns.clustermap(pres, cmap="Blues", cbar_pos=None, linewidths=0.5,
                   linecolor="lightgrey", xticklabels=True, yticklabels=True,
                   figsize=(9, 0.28 * len(pres) + 3))
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=60, ha="right")
g.savefig("results/figures/heatmap_presence_absence.png", dpi=300, bbox_inches="tight")
plt.close("all")

# 2. acquired genes vs mutations per strain
burden = (core.groupby(["strain", "category"]).size().unstack(fill_value=0)
              .reindex(index=strains, columns=["acquired", "mutation"], fill_value=0))
burden = burden.loc[burden.sum(axis=1).sort_values().index]
burden.to_csv("results/burden_per_strain.csv")
burden.plot(kind="barh", stacked=True, figsize=(8, 5), color=["#2b6cb0", "#dd6b20"])
plt.xlabel("Distinct AMR determinants (intrinsic genes excluded)")
plt.title("Acquired genes vs chromosomal mutations per strain")
plt.tight_layout()
plt.savefig("results/figures/burden_per_strain.png", dpi=300)
plt.close("all")

# 3. antibiotic classes per strain
c2 = core.copy()
c2["Class"] = c2["Class"].str.split("/")
c2 = c2.explode("Class")
c2 = c2[c2["Class"] != "MULTIDRUG"]
cls = (c2.groupby(["strain", "Class"]).size().unstack(fill_value=0)
         .reindex(strains, fill_value=0))
cls = cls.loc[(cls > 0).sum(axis=1).sort_values().index]
cls.to_csv("results/class_counts.csv")
cls.plot(kind="barh", stacked=True, figsize=(10, 6), colormap="tab20")
plt.xlabel("Distinct AMR determinants")
plt.title("Resistance determinants by antibiotic class")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
plt.tight_layout()
plt.savefig("results/figures/class_distribution.png", dpi=300)
plt.close("all")

n_classes = (cls > 0).sum(axis=1).rename("n_classes")
n_classes.to_csv("results/classes_per_strain.csv")
print(burden)
print(n_classes)
print("MDR (>=3 classes):", list(n_classes[n_classes >= 3].index))
