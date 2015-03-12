#!/bin/env python3
"""
===================================================
rmdup removes duplicate reads from a sorted,
indexed .sam, .sam.gz, or .bam file.

It uses the python random.choice function to keep
choose the duplicate to keep, thereby reducing the
likelihood of reference bias to almost nothing.

It requires that the pysam package is installed

To run, either run as a script:
    rmdup <infile> <outfile>
or as a module:
    from rmdup import remove_duplicates
    remove_duplicates(infile, outfile, logfile)

The remove_duplicates function requires that the
<infile> and <outfile> are open pysam file handles
and <logfile> is a file name.

Note: The only relevant function in this file is
remove_duplicates(), the rest are all just for
parsing and printing.
==================================================="""
from wasp import mapping

def _get_args():
    """ Command Line Argument Parsing """
    import argparse

    parser = argparse.ArgumentParser(
                 description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input_bam', help="input BAM or SAM file - must be sorted by position and indexed")
    parser.add_argument("output_bam", help="output BAM or SAM file")

    # Optional log file
    parser.add_argument('-l', "--logfile", help="Log file name. Defaults to input file name + .log")

    args = parser.parse_args()

    return(args.input_bam, args.output_bam, args.logfile)

if __name__ == '__main__':
    """Run directly"""
    mapping.remove_duplicates(*_get_args())

