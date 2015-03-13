#!/bin/env python3
"""
===============================================================================
                        WASP - find_intersecting_snps

Search a previously mapped bam or sam file for mapped reads that overlap
single nucleotide polymorphisms (SNPs). For every read that overlaps a SNP, its
genotype is swapped with that of the other allele and the new fastq file is
output to <input_name>.remap.fq.gz.
This file can then be remapped to the genome and filter_remapped_reads can be
used to retrieve properly mapped reads.

Output files:
<input>.sort.bam        - Sorted bamfile of the original input (only produced
                          if sort requested)
<input>.keep.bam        - bamfile with reads that did not intersect SNPs and
                          therefore can be kept without remapping
<input>.to.remap.bam    - bamfile with original reads that need to be remapped
<input>.to.remap.num.gz - matched lines with the input.to.remap.bam that
                          indicate the number of variants of the original read
                          that must be remapped
<input>.remap.fq.gz     - fastq file containing the new variants to remap.
                          Will be .fq1.gz and .fq2.gz if the paired end option
                          is used
===============================================================================
"""
from wasp import mapping
from pysam import sort
import sys, argparse

def _get_args():
    """ Command Line Parsing """
    parser = argparse.ArgumentParser(
                    description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)

    # Files
    parser.add_argument('infile', help='.bam file from the initial mapping process')
    parser.add_argument('snp_dir', help='Directory containing SNPs segregating within the sample in question. Should contain sorted files of SNPs separated by chromosome and named: chr<#>.snps.txt.gz with three columns: position RefAllele AltAllele')

    # Options
    parser.add_argument('-o', dest='outfile_prefix',
                        help='Prefix for outfiles, default is same prefix as infile')
    parser.add_argument('-m', dest='max_window', type=int,
                        default=100000, help="Changes the maximum window to search for SNPs.  The default is 100,000 base pairs.  Reads or read pairs that span more than this distance (usually due to splice junctions) will be thrown out.  Increasing this window allows for longer junctions, but may increase run time and memory requirements.")
    parser.add_argument('-p', dest='is_paired_end', action='store_true',
                        default=False, help='Paired end data (default is single)')
    parser.add_argument('-s', dest='sort', action='store_true',
                        help="Sort the input bam file before using (file must be sorted by coordinate, if it isn't specifying this flag will sort it for you)")

    return(parser)

def generate_names(prefix, paired_end):
    """ Return a tuple of outfile names in the correct order for
        the Bam_scanner class.
        Requires a prefix for the outfile names
        paired_end should be set to True for paired end data, False
        otherwise """
    keep_file_name=pref+".keep.bam"
    remap_name=pref+".to.remap.bam"
    remap_num_name=pref+".to.remap.num.gz"

    # Test for paired end
    if args.is_paired_end:
        fastq_names=[pref+".remap.fq1.gz", pref+".remap.fq2.gz"]
    else:
        fastq_names=[pref+".remap.fq.gz"]

    return(keep_file_name, remap_name, remap_num_name, fastq_names)

def main():
    """ Function for running as a script """
    parser = _get_args()
    args = parser.parse_args()

    # Get infile prefix for outfile creation
    pref = args.outfile_prefix if args.outfile_prefix else '.'.join(infile.split('.')[:-1])

    # If a sort is requested, do it now
    if args.sort:
        new_name = pref + '.sort'
        sort(infile, new_name)
        infile = new_name

    # Generate output file names
    names = generate_names(pref, args.is_paired_end)

    # Create bam scanner object - open all the files
    bam_data = mapping.Bam_scanner(args.is_paired_end, args.max_window, infile, *names, snp_dir=args.snp_dir)
    bam_data.fill_table()

    while not bam_data.end_of_file:
        if args.is_paired_end:
            bam_data.empty_slot_paired()
        else:
            bam_data.empty_slot_single()
        bam_data.fill_table()

# The end
if __name__ == '__main__':
    main()
