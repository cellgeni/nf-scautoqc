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
my_parser.add_argument("--seq", default=None, help="type of sequencing")
args = my_parser.parse_args()

samples = args.samples.split(',')
objects = args.objects.split(',')

ads = []
for obj in objects:
    ads.append(sc.read(obj))

samp_ads = [ad.obs['sampleID'].unique()[0] for ad in ads]
new_order = [samples.index(item) for item in samp_ads]

ads1 = [ads[new_order.index(i)] for i in range(0,len(samples))]

del ads
gc.collect()

pooled_ad = anndata.AnnData.concatenate(
    *ads1, batch_key="sampleID", batch_categories=samples
)

del ads1
gc.collect()

if args.seq == 'single-cell':
    suffixes = ("", "_raw", "_spliced", "_unspliced")
else:
    suffixes = ("", "_raw")

for suffix in suffixes:
    pooled_ad.var[f"n_counts{suffix}"] = pooled_ad.var[
        [f"n_counts{suffix}-{sid}" for sid in samples]
    ].sum(axis=1)

for suffix in suffixes:
    pooled_ad.varm[f"n_counts{suffix}"] = pooled_ad.var[
        [f"n_counts{suffix}-{sid}" for sid in samples]
    ].values

for suffix in suffixes:
    for sid in samples:
        del pooled_ad.var[f"n_counts{suffix}-{sid}"]

for suffix in suffixes:
    pooled_ad.var[f"n_cells{suffix}"] = pooled_ad.var[
        [f"n_cells{suffix}-{sid}" for sid in samples]
    ].sum(axis=1)

for suffix in suffixes:
    pooled_ad.varm[f"n_cells{suffix}"] = pooled_ad.var[
        [f"n_cells{suffix}-{sid}" for sid in samples]
    ].values

for suffix in suffixes:
    for sid in samples:
        del pooled_ad.var[f"n_cells{suffix}-{sid}"]

pooled_ad.write(
    "pooled_postqc.h5ad",
    compression="gzip",
)