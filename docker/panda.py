# /usr/bin/env python

import argparse
from netZooPy.panda import Panda
import pandas as pd
import pickle


def slice_matrix(exprs_filename):
    '''
    Slices matrix into scatter chunks and exports slices to file.
    '''
    exprs_df = pd.read_csv(exprs_filename, index_col = 0, header = 0)
    sample_size = len(exprs_df.columns)
    range_list = []
    for i in range(10):
        start = i
        end = 10*(i + 1)
        if i == 10:
            end += sample_size % 10
        range_list.append((start, end))
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
    slice_matrix(args.exprs)
    run_panda(args, panda_output)


if __name__ == "__main__":
    main()
