# /usr/bin/env python

import argparse
from netZooPy.panda import Panda
import numpy as np
import pandas as pd
import pickle
import sys

# Max number of genes to consider. We take the top N 
# based on their mean expression across all samples
NMAX = 15000

def slice_matrix(exprs_filename, output_filename, num_scatters):
    """Slices samples into bins for LIONESS scatter.

    This writes a file which will tell the various WDL shards which samples
    they should work on.

    Args:
        exprs_filename (str): file name for input expression matrix
        output_filename (str): output file name for scatter slices
        num_scatters (int): number of scatters to (mostly) evenly divide
        the samples into
    """
    try:
        num_scatters = int(num_scatters)
    except ValueError as e:
        sys.stderr.write("ValueError: %s\n" % str(e))
        sys.exit()
    exprs_df = pd.read_table(exprs_filename, index_col = 0, header = 0, sep="\t")
    sample_size = exprs_df.shape[1]
    range_list = []
    quantile = int(sample_size / num_scatters)
    for i in range(num_scatters):
        start = quantile * i + 1 # convert ranges to 1-based index
        end = quantile * (i + 1)
        if i == (num_scatters-1):
            end = sample_size # avoids division remainder issues
        range_list.append((start, end))

    # Export ranges to file as TSV
    with open(output_filename, 'w') as handle:
        for start, end in range_list:
            handle.write("%i\t%i\n" % (start, end))


def run_panda(args, panda_output):
    """Run PANDA and saves as pickle.

    Args:
        args (argparse.ArgumentParser()): passed args from main()
        panda_output (str): output file name for PANDA TSV
    """    
    # Load the data as a pandas dataframes
    # Required for LIONESS
    exprs_df = pd.read_table(args.exprs, index_col = 0, header = 0, sep = "\t")
    motif_df = pd.read_table(args.motif, header = None, sep = "\t")
    ppi_df = pd.read_table(args.ppi, header = None, sep = "\t")
    # Adding headers for the PANDAs obj to read
    motif_df.columns =['source','target','weight']

    # We first subset the expression dataframe to retain only the top NMAX
    # by mean expression. Otherwise, memory consumption is too much.

    # This covers a very fringe case here where the __mean__ column might 
    # already be in the matrix (yes, VERY fringe). 
    # Just keep adding underscores to create a unique column name 
    # for the row-mean values.
    mean_col_name = '__mean__'
    while mean_col_name in exprs_df.columns:
        mean_col_name = '_' + mean_col_name + '_'
    exprs_df[mean_col_name] = exprs_df.apply(lambda x: x.mean(), axis=1)

    # retain only the top NMAX and drop that mean value column since we're done with it.
    exprs_df = exprs_df.nlargest(NMAX, mean_col_name)
    exprs_df.drop(mean_col_name, axis=1, inplace=True)

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
        remove_missing=False, 
        keep_expression_matrix=True, 
        save_memory=False, 
        modeProcess='legacy'
    )

    # Save PANDA object as pickle
    with open("panda_obj.pkl", "wb") as f:
        pickle.dump(panda_obj, f)


def main():
    """Runs PANDA on input expression matrix and produces scatters.
    """    
    # Parse args
    parser = argparse.ArgumentParser(
        description="Runs PANDA on input count matrix."
    )
    parser.add_argument(
        "-s", "--scatter", metavar="STR", required=True,
        default="sample_scatter_ranges.txt",
        help="output name for scatter ranges"
    )
    parser.add_argument(
        "--num_scatters", metavar="INT", required=True,
        help="Number of scatters to make"
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
    slice_matrix(args.exprs, args.scatter, args.num_scatters)
    run_panda(args, args.out)


if __name__ == "__main__":
    main()
