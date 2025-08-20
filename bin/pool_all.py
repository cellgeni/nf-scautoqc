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
ranges = args.ranges.split(',')

if len(samples) != len(objects):
    objects = [i for i in objects if i.startswith(tuple(samples))]

ads = {}
for obj in objects:
    ad = sc.read(obj)
    sample_id = ad.obs['sampleID'].unique()[0]
    ad.var_names_make_unique()
    if sample_id in samples:
        ads[sample_id] = ad

pooled_ad = sc.AnnData.concatenate(
    *[ads[sample] for sample in samples],
    batch_key="sampleID",
    batch_categories=samples
)

wide_series_by_sample = {}

for r in ranges:
    df = pd.read_csv(r, skipinitialspace=True)

    first_col = df.columns[0]
    sample = first_col.split("|", 1)[1]

    s = (
        df.rename(columns={first_col: "metric"})
          .set_index("metric")
          .stack()
    )
    s.index = pd.MultiIndex.from_tuples(s.index, names=["metric", "bound"])
    wide_series_by_sample[sample] = s

wide = pd.DataFrame(wide_series_by_sample).T
wide.index.name = "sample"

wide = wide.sort_index()

metrics = list(dict.fromkeys(wide.columns.get_level_values(0)))  # preserves first-seen order
col_order = pd.MultiIndex.from_product([metrics, ["low", "high"]], names=["metric", "bound"])
wide = wide.reindex(columns=col_order)

wide.to_csv("qc_thresholds.csv")   # columns like ('n_counts','low'), ('n_counts','high'), ...
  
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