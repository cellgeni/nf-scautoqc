#!/usr/bin/env python
"""Add metadata to pooled object

Usage: add_scrublet_meta.py --obj pooled_obj.h5ad --scr /path/to/scrublet_csv1,/path/to/scrublet_csv2...

Options:
  --obj        path of the pooled h5ad
  --scr        path of scrublet outputs
  --meta       csv of metadata of the samples
"""

import scanpy as sc
import pandas as pd
import sctk as sk
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np
import anndata
import gc
import re
import argparse


my_parser = argparse.ArgumentParser()
my_parser.add_argument("--obj", default=None, help="path of the pooled h5ad")
my_parser.add_argument("--scr", default=None, help="path of scrublet outputs")
my_parser.add_argument("--meta", default=None, help="csv of metadata of the samples", nargs='?')
args = my_parser.parse_args()

pooled_ad = sc.read(args.obj)

pooled_ad0 = anndata.AnnData(
    X=pooled_ad.X, obs=pooled_ad.obs.copy(), var=pooled_ad.var.copy()
)

del pooled_ad
gc.collect()

if not args.meta == None:
    metadata = pd.read_csv(args.meta) # assuming sample IDs are in sampleID column

    obs_df = (
        pooled_ad0.obs.reset_index()
        .merge(
            metadata,
            how="left",
            left_on="sampleID",
            right_on="sampleID",
        )
        .set_index("index")
    )

    if not (obs_df.index == pooled_ad0.obs_names).all():
        print('metadata index and object index are not identical, dying...')
        exit()

    pooled_ad0.obs = obs_df

CC_GENE_LIST = "/lustre/scratch127/cellgen/cellgeni/cakirb/QC_pipeline/feature_list/JP_cycle_genes.list"
IG_GENE_LIST = "/lustre/scratch127/cellgen/cellgeni/cakirb/QC_pipeline/feature_list/ig_genes.list"
TCR_GENE_LIST = "/lustre/scratch127/cellgen/cellgeni/cakirb/QC_pipeline/feature_list/tcr_genes.list"

cc_genes = pd.read_csv(CC_GENE_LIST, header=None, names=["gene"]).gene
ig_genes = pd.read_csv(IG_GENE_LIST, header=None, names=["gene"]).gene
tcr_genes = pd.read_csv(TCR_GENE_LIST, header=None, names=["gene"]).gene

igc_genes = [g for g in ig_genes if re.search(r"^IGH[ADEMG][0-9]*$", g)]
ig_genes = ig_genes[~ig_genes.isin(igc_genes)]


pooled_ad0.var["cc"] = pooled_ad0.var_names.isin(cc_genes)
pooled_ad0.var["ig"] = pooled_ad0.var_names.isin(ig_genes)
pooled_ad0.var["tcr"] = pooled_ad0.var_names.isin(tcr_genes)

sk.clear_colors(pooled_ad0)

sample_passqc_df = pd.concat(
    [
        pooled_ad0.obs.groupby("sampleID").apply(
            lambda df: df.good_qc_cluster_mito80.mean()
        ),
        pooled_ad0.obs.groupby("sampleID").apply(
            lambda df: df.good_qc_cluster_mito50.mean()
        ),
        pooled_ad0.obs.groupby("sampleID").apply(
            lambda df: df.good_qc_cluster_mito20.mean()
        ),
        pooled_ad0.obs.groupby("sampleID").apply(
            lambda df: df.good_qc_cluster_mito80.sum()
        ),
        pooled_ad0.obs.groupby("sampleID").apply(
            lambda df: df.good_qc_cluster_mito50.sum()
        ),
        pooled_ad0.obs.groupby("sampleID").apply(
            lambda df: df.good_qc_cluster_mito20.sum()
        ),
        pooled_ad0.obs.sampleID.value_counts(sort=False),
        pooled_ad0.obs[["sampleID"]].drop_duplicates().set_index("sampleID"),
    ],
    axis=1,
).rename(
    columns={
        0: "passQC_frac80",
        1: "passQC_frac50",
        2: "passQC_frac20",
        3: "passQC_count80",
        4: "passQC_count50",
        5: "passQC_count20",
        "count": "total_cell_count",
    }
)

sample_passqc_df.to_csv('sample_passqc_df.csv')

sk.set_figsize((5, 5))
ax = sn.scatterplot(
    data=sample_passqc_df,
    x="total_cell_count",
    y="passQC_count80",
    hue=sample_passqc_df.index,
    size="passQC_frac80",
    sizes=(10, 40),
)
ax.set_xscale("symlog")
ax.set_yscale("symlog")
ax.axline((1e3, 1e3), (1e4, 1e4), c="grey", linestyle="--")
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), ncol=2);
plt.savefig("scatterplot_sample_passqc.png")


sk.set_figsize((4, 5))
fig, axs = plt.subplots(
    nrows=3, ncols=1, sharex=True, sharey=True, constrained_layout=True
)
for i, name in enumerate(["passQC_count80", "passQC_count50", "passQC_count20"]):
    axs[i].hist(np.log10(1 + sample_passqc_df[name].values), 100)
    axs[i].set_title(name)
plt.savefig("barplots_sample_passqc.png")

samples_to_check1 = sample_passqc_df.index.to_series()[
    sample_passqc_df.passQC_count20 > 10000
].to_list()

samples_to_check2 = sample_passqc_df[
    (sample_passqc_df.passQC_frac80 < 0.1) & (sample_passqc_df.passQC_count80 < 100)
].index.to_list()

samples_to_check3 = sample_passqc_df[
    sample_passqc_df.passQC_count80 == 0
].index.to_list()

print(f"samples that have count20 higher than 10k : {samples_to_check1}")
print(f"samples that have frac80 lower than 0.1 & count80 less than 100 : {samples_to_check2}")
print(f"samples that have zero count80: {samples_to_check3}")

pooled_ad1 = pooled_ad0[
    (pooled_ad0.obs.good_qc_cluster <= 80)
    & ~pooled_ad0.obs.sampleID.isin(samples_to_check2)
].copy()

del pooled_ad0
gc.collect()

scr_dfs = []
scr_csvs = args.scr.split(',')
for file in scr_csvs:
    scr_df = pd.read_csv(file, index_col=0)
    scr_df.index = (scr_df.index.to_series().astype(str) + f"-{scr_df.index.name}").values
    scr_dfs.append(scr_df)

scrublet_results_df = pd.concat(scr_dfs)

obs_df = pooled_ad1.obs.merge(
    scrublet_results_df, how="left", left_index=True, right_index=True
)

if not (obs_df.index == pooled_ad1.obs_names).all():
    print('scrublet index and object index are not identical, dying...')
    exit()

pooled_ad1.obs = obs_df

pooled_ad1.obs["doublet"] = (pooled_ad1.obs.scrublet_score > 0.3) | (
    pooled_ad1.obs.bh_pval < 0.65
)
pooled_ad1.obs["stringent_doublet"] = (pooled_ad1.obs.scrublet_score > 0.3) | (
    pooled_ad1.obs.bh_pval < 0.05
)

pooled_ad1.write(
    "pooled.gene_cellbender.good_qc_cluster_mito80.doublet_flagged.h5ad",
    compression="gzip",
)