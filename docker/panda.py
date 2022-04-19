# /usr/bin/env python

import argparse
from netZooPy.panda import Panda
import numpy as np
import pandas as pd
import pickle
import sys


def slice_matrix(exprs_filename, num_ranges):
    '''
    Slices matrix into scatter chunks and exports slices to file.
    '''
    try:
        num_ranges = int(num_ranges)
    except ValueError as e:
        sys.stderr.write("ValueError: %s\n" % str(e))
        sys.exit()
    exprs_df = pd.read_csv(exprs_filename, index_col = 0, header = 0)
    sample_size = exprs_df.shape[0]
    range_list = []
    quantile = int(sample_size / num_ranges)
    for i in num_ranges:
        start = quantile * i + 1 # convert ranges to 1-based index
        end = quantile * (i + 1)
        if i == num_ranges:
            end = sample_size + 1 # avoids division remainder issues
    range_list.append((start, end))
    # Export ranges to file as TSV
    with open("sample_scatter_ranges.txt", 'w') as handle:
        for start, end in range_list:
            handle.write("%i\t%i\n" % (start, end))


def run_panda(args, panda_output):
    '''
    Runs PANDA object creation and exports output to file and pickle.
    '''
    # Load the data as a pandas dataframes
    # Required for LIONESS
    exprs_df = pd.read_csv(args.exprs, index_col = 0, header = 0)
    motif_df = pd.read_csv(args.motif, header = None, sep = "\t")
    ppi_df = pd.read_csv(args.ppi, header = None, sep = "\t")
    # Adding headers for the PANDAs obj to read
    motif_df.columns =['source','target','weight']
    # Running pandas with default expected paramaters
    # save_memory = False results in outputting the PANDA network in edge format
    # save_memory = True results in a matrix format
    # Edge format results in massive memory usage after PANDA finishes
    # LIONESS requires keep_expression_matrix
    # Default modeProcess blows memory stack, so have to use "legacy".
    # Pass the pandas dataframes directly rather than the PATHs.
    panda_obj = Panda(
        exprs_df,
        motif_df,
        ppi_df,
        save_tmp=True,
        save_memory=True,
        remove_missing=False,
        keep_expression_matrix=False
    )
    panda_obj.save_panda_results(panda_output)
    with open("panda_obj.pkl", "wb") as f:
        pickle.dump(panda_obj, f)


def main():
    # Parse args
    parser = argparse.ArgumentParser(
        desc="Runs PANDA on input count matrix."
    )
    parser.add_argument(
        "--ranges", metavar="INT", required=True,
        help="Number of slice ranges"
    )
    parser.add_argument(
        "--motif", metavar="TSV", required=True,
        help="Motif data"
    )
    parser.add_argument(
        "--ppi", metavar="TSV", required=True,
        help="PPI data"
    )
    parser.add_argument(
        "exprs", metavar="TSV",
        help="Expression count matrix"
    )
    args = parser.parse_args()
    panda_output = "panda_output.mtx"
    # Run PANDA
    slice_matrix(args.exprs, args.ranges)
    run_panda(args, panda_output)


if __name__ == "__main__":
    main()
