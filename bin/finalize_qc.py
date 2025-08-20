#!/usr/bin/env python
"""Finalize QC and annotate pooled single-cell object

Usage: finalize_qc.py --obj pooled_obj.h5ad --scr /path/to/scrublet_csv1,/path/to/scrublet_csv2... [--meta sample_metadata.csv] --qc_mode MODE

Options:
    --obj        Path of the pooled h5ad file.
    --scr        Comma-separated paths of Scrublet output CSVs.
    --meta       CSV file with sample metadata (optional).
    --qc_mode    QC mode (original, combined, or other).
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
import os
import argparse

base_dir = os.getenv('BASE_DIR')

my_parser = argparse.ArgumentParser()
my_parser.add_argument("--obj", default=None, help="path of the pooled h5ad")
my_parser.add_argument("--scr", default=None, help="path of scrublet outputs")
my_parser.add_argument("--meta", default=None, help="csv of metadata of the samples", nargs='?')
my_parser.add_argument("--qc_mode", default=None, help="QC mode")

args = my_parser.parse_args()

pooled_ad = sc.read(args.obj)

pooled_ad0 = anndata.AnnData(
    X=pooled_ad.X, obs=pooled_ad.obs.copy(), var=pooled_ad.var.copy()
)

del pooled_ad
gc.collect()

if args.meta is not None:
    metadata = pd.read_csv(args.meta) # assuming sample IDs are in sampleID column

    for i in metadata.columns:
        if i.startswith('Unnamed'):
            del metadata[i]
            continue
        metadata[i] = metadata[i].astype('object')
        metadata[i] = metadata[i].astype('category')

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

CC_GENE_LIST = f"{base_dir}/genes_list/JP_cycle_genes.list"
IG_GENE_LIST = f"{base_dir}/genes_list/ig_genes.list"
TCR_GENE_LIST = f"{base_dir}/genes_list/tcr_genes.list"

cc_genes = pd.read_csv(CC_GENE_LIST, header=None, names=["gene"]).gene
ig_genes = pd.read_csv(IG_GENE_LIST, header=None, names=["gene"]).gene
tcr_genes = pd.read_csv(TCR_GENE_LIST, header=None, names=["gene"]).gene

igc_genes = [g for g in ig_genes if re.search(r"^IGH[ADEMG][0-9]*$", g)]
ig_genes = ig_genes[~ig_genes.isin(igc_genes)]


pooled_ad0.var["cc"] = pooled_ad0.var_names.isin(cc_genes)
pooled_ad0.var["ig"] = pooled_ad0.var_names.isin(ig_genes)
pooled_ad0.var["tcr"] = pooled_ad0.var_names.isin(tcr_genes)

sk.clear_colors(pooled_ad0)

def compute_sample_passqc_df(obs: pd.DataFrame, mode: str) -> pd.DataFrame:
    grp = obs.groupby("sampleID")
    if mode in ("original", "combined"):
        prefix = "good_qc_cluster_mito" if mode == "original" else "consensus_passed_qc_mito"
        df = pd.DataFrame(
            {
                "passQC_frac80": grp[f"{prefix}80"].mean(),
                "passQC_frac50": grp[f"{prefix}50"].mean(),
                "passQC_frac20": grp[f"{prefix}20"].mean(),
                "passQC_count80": grp[f"{prefix}80"].sum(),
                "passQC_count50": grp[f"{prefix}50"].sum(),
                "passQC_count20": grp[f"{prefix}20"].sum(),
                "total_cell_count": grp.size(),
            }
        )
    else:
        df = pd.DataFrame(
            {
                "passQC_frac": grp["consensus_passed_qc"].mean(),
                "passQC_count": grp["consensus_passed_qc"].sum(),
                "total_cell_count": grp.size(),
            }
        )
    df.index.name = "sampleID"
    return df.sort_index()

sample_passqc_df = compute_sample_passqc_df(pooled_ad0.obs, args.qc_mode)
sample_passqc_df.to_csv("sample_passqc_df.csv")

sk.set_figsize((5, 5))

if args.qc_mode in ('original', 'combined'):
    col_count = "passQC_count80"
    col_count2 = "passQC_count20"
    col_frac = "passQC_frac80"
    col_barplots = ["passQC_count80", "passQC_count50", "passQC_count20"]
else:
    col_count = "passQC_count"
    col_count2 = "passQC_count"
    col_frac = "passQC_frac"
    col_barplots = ["passQC_count"]

ax = sn.scatterplot(
    data=sample_passqc_df,
    x="total_cell_count",
    y=col_count,
    hue=sample_passqc_df.index,
    size=col_frac,
    sizes=(10, 40),
)
ax.set_xscale("symlog")
ax.set_yscale("symlog")
ax.axline((1e3, 1e3), (1e4, 1e4), c="grey", linestyle="--")
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), ncol=2)
plt.savefig("scatterplot_sample_passqc.png")
    
sk.set_figsize((4, 5))
fig, axs = plt.subplots(
    nrows=len(col_barplots), ncols=1, sharex=True, sharey=True, constrained_layout=True
)
if not args.qc_mode in ('original', 'combined'):
    axs = np.atleast_1d(axs)
    
for i, name in enumerate(col_barplots):
    axs[i].hist(np.log10(1 + sample_passqc_df[name].values), 100)
    axs[i].set_title(name)
plt.savefig("barplots_sample_passqc.png")

samples_to_check1 = sample_passqc_df.index.to_series()[
    sample_passqc_df[col_count2] > 10000
].to_list()

samples_to_check2 = sample_passqc_df[
    (sample_passqc_df[col_frac] < 0.1) & (sample_passqc_df[col_count] < 100)
].index.to_list()

samples_to_check3 = sample_passqc_df[
    sample_passqc_df[col_count] == 0
].index.to_list()

print(f"samples that have count20 higher than 10k : {samples_to_check1}")
print(f"samples that have frac80 lower than 0.1 & count80 less than 100 : {samples_to_check2}")
print(f"samples that have zero count80: {samples_to_check3}")

scr_dfs = []
scr_csvs = args.scr.split(',')
for file in scr_csvs:
    scr_df = pd.read_csv(file, index_col=0)
    scr_df.index = (scr_df.index.to_series().astype(str) + f"-{scr_df.index.name}").values
    scr_dfs.append(scr_df)

scrublet_results_df = pd.concat(scr_dfs)

obs_df = pooled_ad0.obs.merge(
    scrublet_results_df, how="left", left_index=True, right_index=True
)

obs_df['scrublet_done'] = obs_df['scrublet_done'].astype('category')

if not (obs_df.index == pooled_ad0.obs_names).all():
    print('scrublet index and object index are not identical, dying...')
    exit()

pooled_ad0.obs = obs_df

pooled_ad0.obs["doublet"] = (pooled_ad0.obs.scrublet_score > 0.3) | (
    pooled_ad0.obs.bh_pval < 0.65
)
pooled_ad0.obs["stringent_doublet"] = (pooled_ad0.obs.scrublet_score > 0.3) | (
    pooled_ad0.obs.bh_pval < 0.05
)

pooled_ad0.write(
    "scautoqc_pooled.h5ad",
    compression="gzip",
)

if args.qc_mode == 'original':
    pooled_ad1 = pooled_ad0[
        (pooled_ad0.obs.good_qc_cluster <= 80)
        & ~pooled_ad0.obs.sampleID.isin(samples_to_check2)
    ].copy()
elif args.qc_mode == 'combined':
    pooled_ad1 = pooled_ad0[
        (pooled_ad0.obs.consensus_passed_qc <= 80)
        & ~pooled_ad0.obs.sampleID.isin(samples_to_check2)
    ].copy()
else:
    pooled_ad1 = pooled_ad0[
        (pooled_ad0.obs.consensus_passed_qc == True)
        & ~pooled_ad0.obs.sampleID.isin(samples_to_check2)
    ].copy()

del pooled_ad0
gc.collect()

pooled_ad1.write(
    "scautoqc_pooled_filtered.h5ad",
    compression="gzip",
)