#!/usr/bin/env python
"""Pool post-qc objects

Usage: pool.py --samples sample1,sample2... --objects sample1,sample2...

Options:
  --samples           samples separated by comma
  --objects           object paths separated by comma
"""


import scanpy as sc
import pandas as pd
import anndata
import gc
import argparse

my_parser = argparse.ArgumentParser()
my_parser.add_argument("--samples", default=None, help="samples separated by comma")
my_parser.add_argument("--objects", default=None, help="object paths separated by comma")
args = my_parser.parse_args()

samples = args.samples.split(',')
objects = args.objects.split(',')

ads = []
for sid, obj in zip(samples, objects):
    print(sid)
    ads.append(sc.read(obj))
    ads[-1].obs['sampleID'] = sid

pooled_ad = anndata.AnnData.concatenate(
    *ads, batch_key="sampleID", batch_categories=samples
)

del ads
gc.collect()

for suffix in ("", "_raw", "_spliced", "_unspliced"):
    pooled_ad.var[f"n_counts{suffix}"] = pooled_ad.var[
        [f"n_counts{suffix}-{sid}" for sid in samples]
    ].sum(axis=1)

for suffix in ("", "_raw", "_spliced", "_unspliced"):
    pooled_ad.varm[f"n_counts{suffix}"] = pooled_ad.var[
        [f"n_counts{suffix}-{sid}" for sid in samples]
    ].values

for suffix in ("", "_raw", "_spliced", "_unspliced"):
    for sid in samples:
        del pooled_ad.var[f"n_counts{suffix}-{sid}"]

for suffix in ("", "_raw", "_spliced", "_unspliced"):
    pooled_ad.var[f"n_cells{suffix}"] = pooled_ad.var[
        [f"n_cells{suffix}-{sid}" for sid in samples]
    ].sum(axis=1)

for suffix in ("", "_raw", "_spliced", "_unspliced"):
    pooled_ad.varm[f"n_cells{suffix}"] = pooled_ad.var[
        [f"n_cells{suffix}-{sid}" for sid in samples]
    ].values

for suffix in ("", "_raw", "_spliced", "_unspliced"):
    for sid in samples:
        del pooled_ad.var[f"n_cells{suffix}-{sid}"]

pooled_ad.write(
    "pooled.gene_velo_cellbender.post_qc.h5ad",
    compression="gzip",
)