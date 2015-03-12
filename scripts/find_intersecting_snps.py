#!/bin/env python3
"""
===============================================================================
                            find_intersecting_snps

Search a previously mapped bam or sam file for mapped reads that overlap
single nucleotide polymorphisms (SNPs). For every read that overlaps a SNP, its
genotype is swapped with that of the other allele and the new fastq file is
output to <input_name>.remap.fq.gz.
This file can then be remapped to the genome and filter_remapped_reads can be
used to retrieve properly mapped reads.

Output files:
<input>.sort.bam        - Sorted bamfile of the original input
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
import sys, argparse

def main():
    parser = argparse.ArgumentParser(
                 description=__doc__,
                 formatter_class=argparse.RawDescriptionHelpFormatter)

    # Arguments
    parser.add_argument("-p", action='store_true', dest='is_paired_end',
                        default=False, help="Paired end data (default is single)")
    parser.add_argument("-m", action='store', dest='max_window', type=int,
                        default=100000, help="Changes the maximum window to search for SNPs.  The default is 100,000 base pairs.  Reads or read pairs that span more than this distance (usually due to splice junctions) will be thrown out.  Increasing this window allows for longer junctions, but may increase run time and memory requirements.")
    parser.add_argument("infile", help=".bam file from the initial mapping process")
    parser.add_argument("snp_dir", help="Directory containing SNPs segregating within the sample in question. Should contain sorted files of SNPs separated by chromosome and named: chr<#>.snps.txt.gz with three columns: position RefAllele AltAllele")

    args=parser.parse_args()
    infile=args.infile
    snp_dir=args.snp_dir
    name_split=infile.split(".")

    if len(name_split)>1:
        pref=".".join(name_split[:-1])
    else:
        pref=name_split[0]

    pysam.sort(infile,pref+".sort")

    # Generate output file names
    sort_file_name=pref+".sort.bam"
    keep_file_name=pref+".keep.bam"
    remap_name=pref+".to.remap.bam"
    remap_num_name=pref+".to.remap.num.gz"

    # Test for paired end
    if args.is_paired_end:
        fastq_names=[pref+".remap.fq1.gz",pref+".remap.fq2.gz"]
    else:
        fastq_names=[pref+".remap.fq.gz"]

    bam_data=mapping.bam_scanner(args.is_paired_end,args.max_window,sort_file_name,keep_file_name,remap_name,remap_num_name,fastq_names,snp_dir)
    bam_data.fill_table()
    #i=0
    while not bam_data.end_of_file:
        #i+=1
        #if i>50000:
            #sys.stderr.write(str(asizeof.asizeof(bam_data))+"\t"+str(asizeof.asizeof(bam_data.snp_table))+"\t"+str(asizeof.asizeof(bam_data.read_table))+"\t"+str(bam_data.num_reads)+"\t"+str(bam_data.num_snps)+"\n")
            #sys.stderr.write(str(asizeof.asizeof(bam_data))+"\t"+str(bam_data.num_reads)+"\t"+str(bam_data.num_snps)+"\t"+str(len(bam_data.indel_dict))+"\n")
            #i=0
        if args.is_paired_end:
            bam_data.empty_slot_paired()
        else:
            bam_data.empty_slot_single()
        bam_data.fill_table()

    sys.stderr.write("Finished!\n")

if __name__ == '__main__':
    """Run directly"""
    main()
