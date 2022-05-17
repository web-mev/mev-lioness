#! /usr/bin/env python

import argparse
import pandas as pd

def get_num_samples(fname):
    """Parses TSV for number of columns

    Args:
        fname (str): filename of input TSV

    Returns:
        int: number of columns in TSV
    """    
    df = pd.read_csv(fname, index_col=0, header=0, sep="\t")
    samples = df.columns.values.tolist()
    return len(samples)


def determine_num_ranges(N, max_n):
    """Calculates the number of ranges given max_n from N.

    Args:
        N (int): total N of samples / observations
        max_n (int): max n for a single range

    Returns:
        int: total number of windows
    """    
    return int(N / max_n)


def main():
    """Determines the number of WDL scatters from observations in TSV.
    """    
    parser = argparse.ArgumentParser(
        description="Determines the number of WDL scatters from TSV."
    )
    parser.add_argument(
        "-m", "--max", metavar="INT", default="50",
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
