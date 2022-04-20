#! /usr/bin/env python

import argparse
import pandas as pd
import sys


def order_scatter(list_npys):
    # Order numpy matrices
    ordered_dict = dict(enumerate(list_npys))
    return ordered_dict


def main():
    '''Merges scattered LIONESS matrices into a single target score matrix.'''
    # Parse args
    parser = argparse.ArgumentParser(
        desc="Merges LIONESS matrices into a single target score matrix."
    )
    parser.add_argument(
        "lioness", metavar="NPY", nargs="+",
        help="LIONESS numpy file(s)"
    )
    args = parser.parse_args


if __name__ == "__main__":
    main()
