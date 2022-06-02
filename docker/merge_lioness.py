#! /usr/bin/env python

import argparse
import pandas as pd
import sys


def merge_lioness_scatter(input_files, genes_output, tfs_output):
    """Merges lioness slices into complete matrix and exports to file.

    Args:
        input_files (list): list of input file names
        genes_output (str): file name for gene target score matrix
        tfs_output (str): file name for tf target score matrix
    """    
    # this uses a multiIndex which combines the first two columns
    # containing the transcription factor and gene
    df = pd.read_table(input_files[0], sep="\t", index_col=[0,1])

    # perform the target score sums on both gene and tf.
    # This splits the `df` dataframe so that the joins below
    # are far smaller.
    gene_df = df.groupby(level='gene').sum()
    tf_df = df.groupby(level='tf').sum()
    for fname in input_files[1:]:
        # read and perform the target score summations ahead of time
        other = pd.read_table(fname, sep="\t", index_col=[0,1])
        other_gene_df = other.groupby(level='gene').sum()
        other_tf_df = other.groupby(level='tf').sum()

        # now merge into the 'master' dataframes 
        gene_df = gene_df.join(other_gene_df)
        tf_df = tf_df.join(other_tf_df)

    gene_df.to_csv(
        genes_output,
        sep='\t',
        float_format='%.3f'
    )
    tf_df.to_csv(
        tfs_output,
        sep='\t',
        float_format='%.3f'
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
    merge_lioness_scatter(args.lioness, args.gene, args.tf)


if __name__ == "__main__":
    main()
