#!/usr/bin/env python
"""Add metadata and scrublet scores to pooled single-cell object

Usage: finalize_qc_basic.py --obj pooled_obj.h5ad --scr /path/to/scrublet_csv1,/path/to/scrublet_csv2... [--meta sample_metadata.csv]

Options:
    --obj        Path of the pooled h5ad file.
    --scr        Comma-separated paths of Scrublet output CSVs.
    --meta       CSV file with sample metadata (optional).
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
args = my_parser.parse_args()

pooled_ad = sc.read(args.obj)

pooled_ad0 = anndata.AnnData(
    X=pooled_ad.X, obs=pooled_ad.obs.copy(), var=pooled_ad.var.copy()
)

del pooled_ad
gc.collect()

if not args.meta == None:
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
    "scautoqc_pooled_doubletflagged_metaadded_basic.h5ad",
    compression="gzip",
)