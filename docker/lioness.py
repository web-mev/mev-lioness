# /usr/bin/env python

import argparse
from netZooPy.panda import Panda
from netZooPy.lioness import Lioness
from netZooPy.lioness.analyze_lioness import AnalyzeLioness
import pandas as pd
import pickle
import sys


def run_lioness(panda_obj, start, end):
    '''Runs LIONESS with PANDA object.'''
    # start: start - 1 in python
    # end: end
    # LIONESS start and end is inclusive
    # i.e. translated into R 1-based style
    lioness_obj = Lioness(
        panda_obj, 
        save_dir='../data',
        start = start,
        end = end
    )
    # Save to binary numpy format
    out_filename = "lioness_matrix.npy"
    np.save(
        out_filename,
        np.array(lioness_obj.total_lioness_network),
        allow_pickle=False
    )
    return lioness_obj


def load_panda_obj(panda_filename):
    '''Loads pickle and returns object.'''
    with open(panda_filename, 'r') as f:
        panda_obj = pickle.load(f)
        return panda_obj


def parse_slices(tsv_filename, line_num):
    '''Parses TSV for line number denoting slice ranges.'''
    slice_ranges = []
    with open(tsv_filename, 'r') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            slice_ranges.append(arow)
    return slice_ranges[line_num]


def main():
    # Parse args
    parser = argparse.ArgumentParser(
        desc="Runs LIONESS on input PANDA pickle"
    )
    parser.add_argument(
        "--slices", metavar="TSV", required=True,
        help="TSV of scatter slices"
    )
    parser.add_argument(
        "--line", metavar="INT", required=True,
        help="Line number for slice ranges"
    )
    parser.add_argument(
        "panda", metavar="PICKLE",
        help="PANDA pickle object"
    )
    args = parser.parse_args()
    # Load slicing ranges
    start, end = parse_slices(args.slices, args.line)
    try:
        start, end = int(start), int(end)
    except ValueError as e:
        sys.stderr.write("ValueError: %s\n" % str(e))
        sys.exit(2)
    # Print slicing range to file
    # This is to track sample names
    with open("lioness_slice.txt", "r") as handle:
        handle.write("%i\t%i\n" % (start, end))
    # Load PANDA object
    panda_obj = load_panda_obj(args.panda)
    # Run LIONESS
    run_lioness(panda_obj, start, end)


if __name__ == "__main__":
    main()
