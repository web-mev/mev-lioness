#! /usr/bin/env python

import argparse
import pandas as pd
import sys

def get_num_samples(fname):
    """Parses TSV for number of columns

    Args:
        fname (str): filename of input TSV

    Returns:
        int: number of columns in TSV
    """    
    df = pd.read_table(fname, index_col=0, header=0, sep="\t")
    return df.shape[1]


def determine_num_ranges(N, max_n):
    """Calculates the number of scatters given max_n from N.

    To increase the speed of the LIONESS analysis, we use split the expression
    matrix and operate on smaller subsets (e.g. like a map-reduce).
    To limit the number of samples in a given shard, we specify `max_n`.
    Given that `max_n`, we get the number of shards we will need to make

    Args:
        N (int): total N of samples / observations
        max_n (int): maximum number of samples for a single shard

    Returns:
        int: total number of windows
    """    
    d = int(N / max_n)

    # edge case where max_n > N.
    # We need at least one shard 
    if d == 0:
        return 1
    return d


def main():
    """Determines the number of WDL scatters from observations in TSV.
    """    
    parser = argparse.ArgumentParser(
        description="Determines the number of WDL scatters from TSV."
    )
    parser.add_argument(
        "-m", "--max", metavar="INT", default=50, type=int,
        help="max number of samples per scatter [default: 50]"
    )
    parser.add_argument(
        "tsv", metavar="TSV", help="TSV expression matrix"
    )
    args = parser.parse_args()
    n_ranges = determine_num_ranges(
        get_num_samples(args.tsv),
        args.max
    )
    sys.stdout.write(str(n_ranges))


if __name__ == "__main__":
    main()
