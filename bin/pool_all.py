#!/usr/bin/env python
"""Pool post-qc objects

Usage: pool.py --samples sample1,sample2... --objects sample1,sample2...

Options:
  --samples           samples separated by comma
  --objects           object paths separated by comma
"""


import scanpy as sc
import gc
import argparse

my_parser = argparse.ArgumentParser()
my_parser.add_argument("--samples", default=None, help="samples separated by comma")
my_parser.add_argument("--objects", default=None, help="object paths separated by comma")
my_parser.add_argument("--ss_out", default=None, help="type of sequencing")
args = my_parser.parse_args()

samples = sorted(args.samples.split(','))
objects = args.objects.split(',')

if len(samples) != len(objects):
    objects = [i for i in objects if i.startswith(tuple(samples))]

ads = {}
for obj in objects:
    ad = sc.read(obj)
    sample_id = ad.obs['sampleID'].unique()[0]
    if sample_id in samples:
        ads[sample_id] = ad

pooled_ad = sc.AnnData.concatenate(
    *[ads[sample] for sample in samples],
    batch_key="sampleID",
    batch_categories=samples
)

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
    "scautoqc_pooled.h5ad",
    compression="gzip",
)