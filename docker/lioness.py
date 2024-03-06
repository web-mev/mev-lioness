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
    results.to_csv(output_filename, sep="\t", index=False)


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


def main():
    """Runs LIONESS on input PANDA pickles and given slice range.
    """    
    # Parse args
    parser = argparse.ArgumentParser(
        description="Runs LIONESS on input PANDA pickle and given slice range."
    )
    parser.add_argument(
        "--start", required=True, type=int,
        help="Start index for slice"
    )
    parser.add_argument(
        "--end", required=True, type=int,
        help="End for slice"
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
    
    # Get sample names from slice indices
    samplenames = get_sample_names(args.exprs, args.start, args.end)

    # Load PANDA object
    panda_obj = load_panda_obj(args.panda)

    # Run LIONESS
    run_lioness(panda_obj, 
        args.start, 
        args.end, 
        samplenames,
        args.output, 
        args.save_dir    
    )


if __name__ == "__main__":
    main()
