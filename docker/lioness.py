# /usr/bin/env python

import argparse
from netZooPy.panda import Panda
from netZooPy.lioness import Lioness
from netZooPy.lioness.analyze_lioness import AnalyzeLioness
import pandas as pd
import pickle
import sys


def get_sample_names(exprs_fname, start, end):
    """Loads PANDA input expression matrix to extract sample names from TSV.

    Args:
        exprs_fname (str): PANDA input expression matrix TSV
        start (int): start index (1-based inclusive)
        end (end): end index (1-based inclusive)

    Returns:
        list: strings for sample names
    """    
    df = pd.read_csv(exprs_fname, index_col=0, header=0)
    return list(df.columns[start - 1 : end])


def run_lioness(panda_obj, start, end, samplenames):
    """Runs LIONESS on slice range to output TSV lioness_scatter_output.tsv.

    Args:
        panda_obj (Panda): PANDA object
        start (int): start index for slice (1-based inclusive)
        end (int): end index for slice (1-based inclusive)
        samplenames (list): list of sample names
    """    
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
    #out_filename = "lioness_matrix.npy"
    #np.save(
    #    out_filename,
    #    np.array(lioness_obj.total_lioness_network),
    #    allow_pickle=False
    #)
    results = lioness_obj.export_lioness_results
    results.columns = ["tf", "gene"] + samplenames
    results.to_csv("lioness_scatter_output.tsv", sep="\t")


def load_panda_obj(panda_filename):
    """Loads PANDA pickle to return PANDA object.

    Args:
        panda_filename (str): pickle filename

    Returns:
        Panda: PANDA object
    """    
    with open(panda_filename, 'rb') as f:
        panda_obj = pickle.load(f)
        return panda_obj


def parse_slices(tsv_filename, line_num):
    """Parses TSV for specified line to return slice indices.

    Args:
        tsv_filename (str): TSV file name
        line_num (int): line number in file to extract

    Returns:
        list: [start, end] as 1-based inclusive range
    """    
    slice_ranges = []
    with open(tsv_filename, 'r') as handle:
        for line in handle:
            arow = line.strip('\n').split('\t')
            slice_ranges.append(arow)
    return slice_ranges[line_num]


def main():
    """Runs LIONESS on input PANDA pickles and given slice range.
    """    
    # Parse args
    parser = argparse.ArgumentParser(
        desc="Runs LIONESS on input PANDA pickle and given slice range."
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
        "--exprs", metavar="TSV", required=True,
        help="Input expression matrix to PANDA"
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
    # Get sample names from slice indices
    samplenames = get_sample_names(args.exprs)
    # Load PANDA object
    panda_obj = load_panda_obj(args.panda)
    # Run LIONESS
    run_lioness(panda_obj, start, end)


if __name__ == "__main__":
    main()
