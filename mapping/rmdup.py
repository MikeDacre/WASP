#!/bin/env python
"""===================================================
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
    remove_duplicates(infile, outfile)

The remove_duplicates function requires that the
<infile> and <outfile> are open pysam file handles
==================================================="""
from random import choice
import pysam, sys, os

def remove_duplicates(infile, outfile):
    """ Randomly pick one duplicate only from a sorted, indexed
        sam/bam file.
        infile: a pysam.Samfile open file handle
        outfile: a pysam.Samfile open file handle """

    chr=""
    pos=0
    linelistplus=[]
    linelistminus=[]
    for line in infile:
        if line.rname!=chr or line.pos!=pos:
            # When we reach the end of the duplicates,
            # randomly pick one from the plus strand
            # and one from the minus strand using choice()
            if(len(linelistplus)>0):
                outfile.write(choice(linelistplus))
            if(len(linelistminus)>0):
                outfile.write(choice(linelistminus))

            # Set chr and pos to new coordinates
            chr=line.rname
            pos=line.pos

            # Reset duplicate lists
            linelistplus=[]
            linelistminus=[]

        if(line.flag==0):
            # Add plus strand to the duplicate list
            linelistplus.append(line)
        if(line.flag==16):
            # Add minus strand to the duplicate list
            linelistminus.append(line)
    infile.close()
    outfile.close()

def _printerr(err):
    """ Print a message in red to STDERR """
    print('\033[91m' + err + '\x1b[0m', file=sys.stderr)

def _get_args():
    """ Command Line Argument Parsing """
    import argparse, sys

    parser = argparse.ArgumentParser(
                 description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input_bam', help="input BAM or SAM file - must be sorted by position and indexed")
    parser.add_argument("output_bam", help="output BAM or SAM file")

    options = parser.parse_args()

    if options.input_bam.endswith(".sam") or options.input_bam.endswith("sam.gz"):
        infile = pysam.Samfile(options.input_bam, "r")
    elif options.input_bam.endswith(".bam"):
        # assume binary BAM file
        infile = pysam.Samfile(options.input_bam, "rb")
    else:
        parser.print_help()
        _printerr("Input file must be a .sam, .sam.gz, or .bam file")
        sys.exit(1)

    if options.output_bam.endswith(".sam"):
        # output in text SAM format
        outfile = pysam.Samfile(options.output_bam, "w", template=infile)
    elif options.output_bam.endswith(".bam"):
        # output in binary compressed BAM format
        outfile = pysam.Samfile(options.output_bam, "wb", template=infile)
    else:
        parser.print_help()
        _printerr("Output file must be a .sam, .sam.gz, or .bam file")
        sys.exit(1)

    return(infile, outfile)

# Main function for direct running
def main():
    """Run directly"""
    # Get commandline arguments
    infile, outfile = _get_args()
    remove_duplicates(infile, outfile)

# The end
if __name__ == '__main__':
    main()

