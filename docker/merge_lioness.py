#! /usr/bin/env python

import argparse
import pandas as pd
import sys


def merge_lioness_scatter(input_files, full, genes, tfs):
    """Merges lioness slices into complete matrix and exports to file.

    Args:
        input_files (list): list of input file names
        full (str): file name for unrolled matrix
        genes (str): file name for gene target score matrix
        tfs (str): file name for tf target score matrix
    """    
    df = pd.read_csv(input_files[0], sep="\t", index_col=0)
    for fname in input_files[1:]:
        other = pd.read_csv(fname, sep="\t", index_col=0)
        df = df.join(other.iloc[:, 2:])
    # Drop gene and tf columns, set index to a combine unique gene<->tf index
    # export to TSV
    df.iloc[:, 2:].set_index(
        df['gene'] + "<->" + df['tf']
    ).to_csv(full, sep="\t")
    num_genes = len(df['gene'].drop_duplicates())
    num_tfs = len(df['tf'].drop_duplicates())
    num_samples = len(df.columns[2:])
    arr = np.array(df.iloc[:, 2:]).reshape(num_genes, num_tfs, num_samples)
    # Sum by gene all tfs, export to file
    pd.DataFrame(
        arr.sum(axis=1),
        index=list(df['gene'].drop_duplicates()),
        columns=df.columns[2:]
    ).to_csv(genes, sep="\t")
    pd.DataFrame(
        arr.sum(axis=0),
        index=list(df['tf'].drop_duplicates()),
        columns=df.columns[2:]
    ).to_csv(tfs, sep="\t")


def main():
    '''Merges scattered LIONESS matrices into a single target score matrix.'''
    # Parse args
    parser = argparse.ArgumentParser(
        desc="Merges LIONESS matrices into a single target score matrix."
    )
    parser.add_argument(
        "--full", metavar="TSV",required=True, 
        help="Output file name for unrolled matrix"
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
