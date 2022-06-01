import argparse
from netZooPy.lioness import Lioness
import pandas as pd
import pickle


def get_sample_names(exprs_fname, start, end):
    """Loads PANDA input expression matrix to extract sample names from TSV.

    Args:
        exprs_fname (str): PANDA input expression matrix TSV
        start (int): start index (1-based inclusive)
        end (end): end index (1-based inclusive)

    Returns:
        list: strings for sample names
    """    
    df = pd.read_table(exprs_fname, index_col=0, header=0, sep='\t')
    # note that the start and end are 1-based, so start-1:end picks up
    # the correct indices in python's 0-based indexing.
    return list(df.columns[start - 1 : end])


def run_lioness(panda_obj, 
    start, 
    end, 
    samplenames, 
    output_filename,
    save_dir):
    """Runs LIONESS on slice range to output TSV lioness_scatter_output.tsv.

    Args:
        panda_obj (Panda): PANDA object
        start (int): start index for slice (1-based inclusive)
        end (int): end index for slice (1-based inclusive)
        samplenames (list): list of sample names
        output_filename (str): name of the output file.
    """    
    # start: start - 1 in python
    # end: end
    # LIONESS start and end is inclusive
    # i.e. translated into R 1-based style
    lioness_obj = Lioness(
        panda_obj, 
        save_dir=save_dir,
        start = start,
        end = end
    )

    results = lioness_obj.export_lioness_results
    results.columns = ["tf", "gene"] + samplenames
    results.to_csv(output_filename, sep="\t")


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
    slice_ranges_df = pd.read_table(
        tsv_filename,
        sep='\t',
        header=None,
        names = ['start', 'end']
    )
    return slice_ranges_df.iloc[line_num]


def main():
    """Runs LIONESS on input PANDA pickles and given slice range.
    """    
    # Parse args
    parser = argparse.ArgumentParser(
        description="Runs LIONESS on input PANDA pickle and given slice range."
    )
    parser.add_argument(
        "--slices", metavar="TSV", required=True,
        help="TSV of scatter slices"
    )
    parser.add_argument(
        "--line", metavar="INT", required=True, type=int,
        help="Line number for slice ranges"
    )
    parser.add_argument(
        "--exprs", metavar="TSV", required=True,
        help="Input expression matrix to PANDA"
    )
    parser.add_argument(
        "--output", metavar="TSV", required=True,
        help="Output filename"
    )
    parser.add_argument(
        "--save_dir", metavar="DIR", required=True,
        help="Directory to save intermediate calcs"
    )
    parser.add_argument(
        "panda", metavar="PICKLE",
        help="PANDA pickle object"
    )
    args = parser.parse_args()
    
    # Load slicing ranges and get the start and end
    # locations for the slice. Note that the numbers
    # are inclusive (e.g. if the end of the range is 
    # 10, then we include that, unlike the range func)
    start, end = parse_slices(args.slices, args.line)
    
    # Get sample names from slice indices
    samplenames = get_sample_names(args.exprs)

    # Load PANDA object
    panda_obj = load_panda_obj(args.panda)

    # Run LIONESS
    run_lioness(panda_obj, 
        start, 
        end, 
        samplenames,
        args.output, 
        args.save_dir    
    )


if __name__ == "__main__":
    main()
