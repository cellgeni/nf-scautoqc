#!/usr/bin/env python
"""Flag doublet by scrublet

Usage: program [options] --samp <sampleid> --input <input>

Options:
  --debug           print debug information
  --profile         print profile information
  --samp            sample id
  --input <input>   input h5ad
"""


import logging
import signal
import sys
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import scanpy as sc
import sctk as sk

import argparse

def main(args):
    logging.debug(args)

    input_h5ad = args.input

    ad = sc.read(input_h5ad)
    if (ad.X.sum(0) > 0).sum() < ad.shape[1] * 0.2:
        exit()

    df = sk.run_scrublet(ad, inplace=False)

    df.index.name = args.samp
    df.to_csv(f"{args.samp}_scrublet.csv")

    return 0


if __name__ == '__main__':
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument("--input", default=None, help="input h5ad")
    my_parser.add_argument("--samp", default=None, help="input h5ad")
    my_parser.add_argument("--filter", default=None, help="only consider cells that have True values in this variable of obs")
    my_parser.add_argument("--groupby", default=None, help="run scrublet for each group specified by this variable of obs")
    my_parser.add_argument("--debug", default=None, help="print debug information")
    my_parser.add_argument("--profile", default=None, help="print profile information")

    args = my_parser.parse_args()

    try:
        if args.debug:
            logLevel = logging.DEBUG
        else:
            logLevel = logging.WARN
        logging.basicConfig(
            level=logLevel,
            format="%(asctime)s; %(levelname)s; %(funcName)s; %(message)s",
            datefmt="%y-%m-%d %H:%M:%S"        )
        if args.profile:
            import cProfile
            cProfile.run("main(args)")
        else:
            main(args)
    except KeyboardInterrupt:
        logging.warning("Interrupted")
        sys.exit(1)
