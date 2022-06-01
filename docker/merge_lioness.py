#! /usr/bin/env python

import argparse
import pandas as pd
import sys


def merge_lioness_scatter(input_files, full, genes_output, tfs_output):
    """Merges lioness slices into complete matrix and exports to file.

    Args:
        input_files (list): list of input file names
        full (str): file name for unrolled matrix
        genes_output (str): file name for gene target score matrix
        tfs_output (str): file name for tf target score matrix
    """    
    # this uses a multiIndex which combines the first two columns
    # containing the transcription factor and gene
    df = pd.read_table(input_files[0], sep="\t", index_col=[0,1])
    for fname in input_files[1:]:
        other = pd.read_table(fname, sep="\t", index_col=[0,1])
        df = df.join(other, how='inner')

    # For each gene, sum across the transcription factors on a 
    # per sample basis to get the targeting score of that gene
    # for that sample:
    df.groupby(level='gene').sum().to_csv(
        genes_output,
        sep='\t'
    )

    # For each TF, sum across the genes on a per sample
    # basis to get the targeting score of that TF
    # for that sample:
    df.groupby(level='tf').sum().to_csv(
        tfs_output,
        sep='\t'
    )


def main():
    '''Merges scattered LIONESS matrices into a single target score matrix.'''
    # Parse args
    parser = argparse.ArgumentParser(
        description="Merges LIONESS matrices into a single target score matrix."
    )
    parser.add_argument(
        "--gene", metavar="TSV", required=True,
        help="Output file name for gene target score matrix"
    )
    parser.add_argument(
        "--tf", metavar="TSV", required=True,
        help="Output file name for tf target score matrix"
    )
    parser.add_argument(
        "--lioness", metavar="TSV", nargs="+", required=True,
        help="unrolled LIONESS TSV(s)"
    )
    args = parser.parse_args()
    merge_lioness_scatter(
        args.lioness, args.full, args.gene, args.tf
    )


if __name__ == "__main__":
    main()
