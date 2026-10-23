# scRNA-seq analysis for BRCA2 deficient prostate cancer pdx
# 2026/8/25, Masahiko Ajiro

import scanpy as sc
import anndata as ad
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sc.logging.print_versions()
PackageVersionscvelo0.3.4session_info1.0.1 (v1.0.1)ComponentInfoPython3.10.20 
DependencyVersiontexttable1.7.0debugpy1.8.21tornado6.5.7setuptools82.0.1python-dateutil2.9.0.post0numba0.65.1pytz2026.2stack_data0.6.3pure_eval0.2.3h5py3.16.0charset-normalizer3.4.7llvmlite0.47.0prompt_toolkit3.0.52asttokens3.0.1six1.17.0jedi0.20.0executing2.2.1decorator5.3.1ipython8.37.0pillow12.2.0natsort8.4.0pickleshare0.7.5cycler0.12.1parso0.8.7igraph0.11.9kiwisolver1.5.0joblib1.5.3patsy1.0.2appnope0.1.4statsmodels0.14.6psutil7.2.2

a = sc.read_10x_h5("h5d/PDX21_filtered.h5")
b = sc.read_10x_h5("h5d/PDX22_filtered.h5")
c = sc.read_10x_h5("h5d/PDX36_23188-19-06-01-01_S28_filtered.h5")

for ad in [a,b,c]:ad.obs_names_make_unique()
for ad in [a,b,c]:ad.var_names_make_unique()
for ad in [a,b,c]:ad.layers["counts"] = ad.X.copy()

for ad in [a,b,c]:sc.pp.normalize_total(ad, target_sum=1e4)
for ad in [a, b, c]:
    sc.pp.calculate_qc_metrics(ad, inplace=True)
for name, x in zip(["drugA","drugB","drugC"], [a, b, c]):
    print(name, x.obs["total_counts"].median())  
for ad in [a,b,c]:sc.pp.log1p(ad)
for ad in [a,b,c]: sc.pp.calculate_qc_metrics(ad, inplace=True)
for x, y in zip(["drugA","drugB","drugC"], [a, b, c]):
    print(x, y.obs["total_counts"].median())
tmp={}
tmp["control"] = a.copy()
tmp["short"] = b.copy()
tmp["long"] = c.copy()
total = ad.concat(tmp, label = "treatment")

sc.tl.pca(total, n_comps = 50, svd_solver = "all")
sc.pp.neighbors(total)
sc.tl.umap(total)
sc.tl.leiden(total, flavor="igraph", n_iterations=-1, resolution=2)
total.var['mt'] = total.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(total, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
t2 = total[total.obs['pct_counts_mt'] < 3, :].copy()

# projection of UMAP and Treatment Conditions
sc.pl.umap(total, size = 1, color = ["leiden", "treatment"])
sc.pl.umap(t2, size = 1, color = ["leiden", "treatment"])

t3 = t2[~t2.obs["leiden"].isin(["25", "27", "30"])].copy()
a = sc.read_10x_h5("h5d/PDX21_filtered.h5")
b = sc.read_10x_h5("h5d/PDX22_filtered.h5")
c = sc.read_10x_h5("h5d/PDX36_23188-19-06-01-01_S28_filtered.h5")

for ad in [a,b,c]:ad.obs_names_make_unique()
for ad in [a,b,c]: ad.var_names_make_unique()
for ad in [a,b,c]: ad.layers["counts"] = ad.X.copy()

a2 = a[a.obs_names.isin(t3.obs_names)].copy()
b2 = b[b.obs_names.isin(t3.obs_names)].copy()
c2 = c[c.obs_names.isin(t3.obs_names)].copy()

for d in [a2,b2,c2]:sc.pp.normalize_total(d, target_sum=1e4)
for d in [a2,b2,c2]:sc.pp.log1p(d)
tmp={}
tmp["control"] = a2.copy()
tmp["short"] = b2.copy()
tmp["long"] = c2.copy()
import anndata as ad
t4 = ad.concat(tmp, label = "treatment")
sc.tl.pca(t4, n_comps = 50, svd_solver = "all")

sc.pp.neighbors(t4)
sc.tl.umap(t4)
sc.tl.leiden(t4, flavor="igraph", n_iterations=-1, resolution=0.5)

my_palette = [
    '#ff99ff', '#cc99ff', '#9999ff', '#99ccff', '#99ffff',
    '#99ff99', '#99ff33', '#ffff00', '#ff9999', '#ff9933'
]
t4.uns["leiden_colors"] = my_palette[:t4.obs["leiden"].nunique()]
my_palette2 = [
    '#00cccc', '#cccc00', '#ff007f'
]
t4.uns["treatment_colors"] = my_palette2[:t4.obs["treatment"].nunique()]

# UMAP and treatment condition projection
sc.pl.umap(t4, size = 1.5, color = "leiden",
           palette=my_palette[:len(t4.obs["leiden"].cat.categories)],
           legend_loc="on data")
sc.pl.umap(t4, size = 1.5, color = ["treatment"])

# analysis for BRCA2 positive cell only
a = sc.read_10x_h5("h5d/PDX21_filtered.h5")
b = sc.read_10x_h5("h5d/PDX22_filtered.h5")
c = sc.read_10x_h5("h5d/PDX36_23188-19-06-01-01_S28_filtered.h5")

for ad in [a,b,c]:ad.obs_names_make_unique()
for ad in [a,b,c]: ad.var_names_make_unique()
for ad in [a,b,c]: ad.layers["counts"] = ad.X.copy()

a3 = a[a[:,"BRCA2"].X > 0].copy()
b3 = b[b[:,"BRCA2"].X > 0].copy()
c3 = c[c[:,"BRCA2"].X > 0].copy()

for d in [a3,b3,c3]:sc.pp.normalize_total(d, target_sum=1e4)
for d in [a3,b3,c3]:sc.pp.log1p(d)
tmp={}
tmp["control"] = a3.copy()
tmp["short"] = b3.copy()
tmp["long"] = c3.copy()
t5 = ad.concat(tmp, label = "treatment")
sc.tl.pca(t5, n_comps = 50, svd_solver = "all")
sc.pp.neighbors(t5)
sc.tl.umap(t5)
sc.tl.leiden(t5, flavor="igraph", n_iterations=-1, resolution=0.5)

# UMAOP projection

sc.pl.umap(t5, size = 5, color = ["leiden"], legend_loc="on data")
sc.pl.umap(t5, size = 2, color = ["treatment"])

# pseudotime violinplot
t5.uns["iroot"] = np.flatnonzero(t5.obs["leiden"] == "0")[0]
sc.tl.dpt(t5)
sc.pl.violin(
    t5,
    keys="dpt_pseudotime",
    groupby="leiden",
    stripplot=False
)

# plot for BRCA2 or POLQ positive cell rates

cluster_map = {
    "0": "L0_2_3_4_5_7_8",
    "2": "L0_2_3_4_5_7_8",
    "3": "L0_2_3_4_5_7_8",
    "4": "L0_2_3_4_5_7_8",
    "5": "L0_2_3_4_5_7_8",
    "7": "L0_2_3_4_5_7_8",
    "8": "L0_2_3_4_5_7_8",
    "1": "L1_6",
    "6": "L1_6",
    "9": "L9_10",
    "10": "L9_10"
}

group_order = [
    "control_total",
    "control_L0_2_3_4_5_7_8",
    "control_L1_6",
    "short_total",
    "short_L0_2_3_4_5_7_8",
    "short_L1_6",
    "long_L9_10"
]

status_order = [
    "BRCA2- POLQ-",
    "BRCA2- POLQ+",
    "BRCA2+ POLQ-",
    "BRCA2+ POLQ+"
]

clusters = list(cluster_map.keys())
adata_sub = t4[
    t4.obs["leiden"].astype(str).isin(clusters)
    & t4.obs["treatment"].isin(["control", "short", "long"])
].copy()

adata_sub.obs["leiden_group"] = (
    adata_sub.obs["leiden"]
    .astype(str)
    .map(cluster_map)
)
df = sc.get.obs_df(
    adata_sub,
    keys=["BRCA2", "POLQ", "treatment", "leiden_group"]
)

df = df[
    (
        df["treatment"].isin(["control", "short"])
        & df["leiden_group"].isin([
            "L0_2_3_4_5_7_8",
            "L1_6"
        ])
    )
    |
    (
        (df["treatment"] == "long")
        & (df["leiden_group"] == "L9_10")
    )
].copy()

df["status"] = np.select(
    [
        (df["BRCA2"] > 0) & (df["POLQ"] == 0),
        (df["BRCA2"] == 0) & (df["POLQ"] > 0),
        (df["BRCA2"] > 0) & (df["POLQ"] > 0),
        (df["BRCA2"] == 0) & (df["POLQ"] == 0)
    ],
    [
        "BRCA2+ POLQ-",
        "BRCA2- POLQ+",
        "BRCA2+ POLQ+",
        "BRCA2- POLQ-"
    ],
    default="NA"
)

df = df[df["status"] != "NA"].copy()
df_grouped = df.copy()
df_grouped["plot_group"] = (
    df_grouped["treatment"].astype(str)
    + "_"
    + df_grouped["leiden_group"].astype(str)
)

df_total = df[
    df["treatment"].isin(["control", "short"])
].copy()

df_total["plot_group"] = (
    df_total["treatment"].astype(str)
    + "_total"
)

plot_df = pd.concat(
    [
        df_grouped,
        df_total
    ],
    axis=0,
    ignore_index=True
)
count_df = pd.crosstab(
    plot_df["plot_group"],
    plot_df["status"]
)
count_df = count_df.reindex(
    index=group_order,
    columns=status_order,
    fill_value=0
)
ratio_df = count_df.div(
    count_df.sum(axis=1).replace(0, np.nan),
    axis=0
).fillna(0)

ax = ratio_df.plot(
    kind="bar",
    stacked=True,
    figsize=(6.5, 4.5),
    width=0.8,
    color={
        "BRCA2+ POLQ-": "#ff00ff",
        "BRCA2- POLQ+": "#99ccff",
        "BRCA2+ POLQ+": "#ff0000",
        "BRCA2- POLQ-": "#a0a0a0"
    },
    edgecolor="black"
)

ax.set_xlabel("")
ax.set_ylabel("Fraction")
ax.set_ylim(0, 1)

ax.legend(
    title="",
    frameon=False,
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)
for i, (_, row) in enumerate(ratio_df.iterrows()):
    bottom = 0
    for status in status_order:
        value = row[status]
        #if value > 0.03:
        if value > 0:
            ax.text(
                i,
                bottom + value / 2,
                f"{value * 100:.1f}%",
                ha="center",
                va="center",
                #fontsize=6,
                fontsize=3,
                color="black"
            )

        bottom += value

plt.xticks(
    ticks=range(len(group_order)),
    labels=group_order,
    rotation=45,
    ha="right"
)

plt.tight_layout()
plt.show()

# gene expression heatmap

sc.tl.rank_genes_groups(
    t4,
    groupby="leiden",
    method="wilcoxon"
)
sc.tl.dendrogram(t4, groupby="leiden")
sc.pl.rank_genes_groups_heatmap(
    t4,
    n_genes=100,
    groupby="leiden",
    swap_axes=True,
    show_gene_labels=False,
    cmap="Reds",
    figsize=(6,4)
)

# heatmap projection for cellular signaling

sv = sc.pl.stacked_violin(
    t4,
    var_names=signatures,
    groupby="leiden",
    standard_scale="var",
    cmap="Reds",
    swap_axes=True,
    figsize=(5.5, 3.5),
    show=False
)
ax = sv["mainplot_ax"]
ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")
plt.setp(
    ax.get_xticklabels(),
    rotation=90,
    ha="left"
)
plt.savefig(
    "stacked_violin_signatures.pdf",
    dpi=450,
    bbox_inches="tight"
)
plt.show()
plt.close()


# heatmap projection by gene signature

signature_genes = {
    "IFNA_signature": ifa_sig, # or IFNG signature : ifg_sig
}

color_map = {
    "control_p": "blue",
    "control_n": "lightblue",
    "short_p": "purple",
    "short_n": "violet",
    "long_p": "red",
    "long_n": "pink"
}

X = t4[:, ifa_sig].to_df()

X[group_col] = test_2.obs[group_col].values  # set group_col

X = X.sort_values(group_col)

row_colors = X[group_col].map(color_map)

X_mat = X.drop(columns=group_col)

g = sns.clustermap(
    X_mat.T,
    row_cluster=True,
    col_cluster=False,
    col_colors=row_colors,
    cmap="afmhot",
    xticklabels=False,
    yticklabels=True,
    figsize=(8, 8)
)

g.ax_heatmap.tick_params(axis="y", labelsize=5)
g.ax_heatmap.set_xlabel("Cells", fontsize=8)
g.ax_heatmap.set_ylabel("Genes", fontsize=8)

plt.show()
