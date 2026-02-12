#!/usr/bin/env python
"""Pool post-qc objects

Usage: pool.py --samples sample1,sample2... --objects sample1,sample2...

Options:
  --samples           samples separated by comma
  --objects           object paths separated by comma
  --ranges            ranges csv paths separated by comma
"""


import scanpy as sc
import pandas as pd
import gc
import argparse

my_parser = argparse.ArgumentParser()
my_parser.add_argument("--samples", default=None, help="samples separated by comma")
my_parser.add_argument("--objects", default=None, help="object paths separated by comma")
my_parser.add_argument("--ranges", default=None, help="ranges csv paths separated by comma")
args = my_parser.parse_args()

samples = sorted(args.samples.split(','))
objects = args.objects.split(',')
if args.ranges != 'subset.py':
    ranges = args.ranges.split(',')
else:
    ranges = None

if len(samples) != len(objects):
    objects = [i for i in objects if i.startswith(tuple(samples))]

ads = {}
for obj in objects:
    ad = sc.read(obj)
    sample_id = ad.obs['sampleID'].unique()[0]
    if len(ad.var_names[0].split('_')) > 1:
        ad.var_names = [f"{j}" for i, j in ad.var_names.str.split('_').str[0:2]]
    ad.var_names_make_unique()
    if sample_id in samples:
        ads[sample_id] = ad

pooled_ad = sc.AnnData.concatenate(
    *[ads[sample] for sample in samples],
    batch_key="sampleID",
    batch_categories=samples
)

wide_series_by_sample = {}

if ranges:
    long_frames = []
    for r in ranges:
        if not r.endswith("_qc_thresholds.csv"):
            continue
        try:
            long_df = pd.read_csv(r, skipinitialspace=True)
            long_frames.append(long_df)
        except FileNotFoundError:
            continue
    if long_frames:
        qc_long = pd.concat(long_frames, ignore_index=True)
        sort_cols = [c for c in ("sample", "sampleID") if c in qc_long.columns]
        if sort_cols:
            qc_long = qc_long.sort_values(by=sort_cols, kind="mergesort")
        qc_long.to_csv("qc_thresholds.csv", index=False)
  
ads.clear()
gc.collect()

suffixes = tuple([""] + [f"_{i}" for i in pooled_ad.layers.keys() if i != 'ambiguous'])

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
    "scautoqc_pooled0.h5ad",
    compression="gzip",
)