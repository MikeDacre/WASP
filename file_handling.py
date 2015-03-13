#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8 tabstop=4 expandtab shiftwidth=4 softtabstop=4
"""
================================================================================
                            WASP - file_handing

Simple functions for handing files, IO, and logging

 Last modified: 2015-03-12 20:15

================================================================================
"""
import pysam, gzip, sys, os

def open_infile(infile):
    """ Take an input sam, sam.gz, or bam file and test it.
        Return pysam file handle """
    try:
        if infile.endswith(".sam") or infile.endswith("sam.gz"):
            pyfile = pysam.Samfile(infile, "r")
        elif infile.endswith(".bam"):
            # assume binary BAM file
            pyfile = pysam.Samfile(infile, "rb")
        else:
            raise ValueError("File does not end in .sam, .sam.gz, or .bam")
    except ValueError as error:
        parser.print_help()
        _printerr("\nInput file error. Cannot open.\n")
        print("Error message:", file=sys.stderr)
        print(error, file=sys.stderr)
        wasp.sys.exit(1)
    except OSError:
        _printerr("Input file {} not found".format(infile))
        wasp.sys.exit(1)
    return(pyfile)

def open_outfile(outfile, template=''):
    """ Take an input file location and open it as a pysam handle
        Optionally you can provide an open pysam input file bam,
        which can then be used as the template for the output file
        Exit with error code 1 on fail """
    if outfile.endswith(".sam"):
        # output in text SAM format
        if template:
            pyfile = pysam.Samfile(outfile, "w", template=template)
        else:
            pyfile = pysam.Samfile(outfile, "w")
    elif outfile.endswith(".bam"):
        # output in binary compressed BAM format
        if template:
            pyfile = pysam.Samfile(outfile, "wb", template=template)
        else:
            pyfile = pysam.Samfile(outfile, "wb")
    else:
        parser.print_help()
        file_handing.printerr("Output file must be a .sam or .bam file")
        wasp.sys.exit(1)

    return(pyfile)

def parse_snpfile(input, verbose=False):
    """ Create a dictionary of SNP positions, base 0, from
        an input. (MUST BE A LIST, one item is fine)

        Possible input file(s):
            .bed, .bed.gz - Assumed that first coordinate is base0
            .vcf, .vcf.gz - Assumed that first coordinate is base1
            chr<##>.snps.txt, chr<##>.snps.txt.gz - This file has the format:
                position RefAllele AltAllele
            Assumed that position is base0

        You can also provide a list of files, or a directory containing
        files (only provide one directory, others ignored)

        Note: file name matters. It is used for parsing

        Returns a dictionary formatted like this:
            {chr: {position: (ref, alt)}}

        POSITION IS BASE0 """
    #################################
    #    File parsing functions     #
    #################################
    def parse_file(file, file_type, gzipped):
        """ Choose a parsing function based on file_type,
            open the file, and pass the filehandle to the correct
            function. Returns the snp dictionary """
        # Open the file
        file = gzip.open(file, 'rb') if gzipped else open(file, 'rb')

        if file_type == 'vcf':
            snps = vcf_parse(file)
        elif file_type == 'bed':
            snps = bed_parse(file)
        elif file_type == 'snpfile':
            snps = snpfile_parse(file)
        else:
            raise Exception(file.name + ' is not of the correct type')

        file.close()
        return(snps)

    def vcf_parse(file):
        """ Return dictionary from open vcf 'rb' filehandle
            Note: I convert to base0 here """
        count = 0
        skipped_snps = 0
        snps= {}
        for i in file:
            i = i.decode().rstrip()
            if i.startswith('#'):
                continue
            count += 1
            if count > 10:
                break
            f = i.split('\t')
            if len(f[3]) > 1 or len(f[4]) > 1:
                skipped_snps += 1
                continue

            if f[0] in snps:
                snps[str(f[0])][int(f[1]) - 1] = (f[3], f[4])
            else:
                snps[str(f[0])] = {int(f[1]) - 1: (f[3], f[4])}

        # Print skipped SNPs
        if skipped_snps:
            print('File {} has {} skipped snps'.format(file.name, skipped_snps), file=sys.stderr)
        return(snps)

    def bed_parse(file):
        """ Return dictionary from open bed 'rb' filehandle """
        skipped_snps = 0

        for i in file:
            i = i.decode().rstrip()
            if i.startswith('#'):
                continue
            f = i.split('\t')
            a = f[3].split('/')
            l = len(a)
            if l == 1:
                a = f[3].split('|')
            if l == 1:
                printerr("Bed file alleles need to be separated by either a / or a |. Please try again")
                sys.exit(2)
            if l > 2:
                skipped_snps += 1
                continue
            if len(a[0]) > 1 or len(a[1]) > 1:
                skipped_snps += 1
                continue

            if f[0] in snps:
                snps[str(f[0])][int(f[1])] = (a[0], a[1])
            else:
                snps[str(f[0])] = {int(f[1]): (a[0], a[1])}

        # Print skipped SNPs
        if skipped_snps:
            print('File {} has {} skipped snps'.format(file.name, skipped_snps), file=sys.stderr)
        return(snps)

    def snpfile_parse(file):
        """ Return dictionary from open snp text file 'rb' filehandle """
        from re import findall
        # Test that chromosome name is in the file Name
        names = file.name.split('.')
        if len(names) < 3:
            printerr("{} has an improperly formatted filename".format(file.name))
            printerr("Required format: chr<#>.snps.txt[.gz]")
            raise TypeError("Improper name")

        snps = {}
        skipped_snps = 0
        chr = findall(r'[0-9]+$', names[0])[0]
        for i in file:
            i = i.decode().rstrip()
            if i.startswith('#'):
                continue
            f = i.split('\t')
            if len(f[1]) > 1 or len(f[2]) > 1:
                skipped_snps += 1
                continue
            if names[0] in snps:
                snps[str(chr)][int(f[0])] = (f[1], f[2])
            else:
                snps[str(chr)] = {int(f[0]): (f[1], f[2])}

        # Print skipped SNPs
        if skipped_snps:
            print('File {} has {} skipped snps'.format(file.name, skipped_snps), file=sys.stderr)
        return(snps)

    #################################
    #         Name parsing          #
    #################################
    def name_parse(file):
        """ Check a file to see what type it is from the filename """
        f = file.split('.')
        gzipped = True if f[-1] == 'gz' else False
        suffix = f[-2] if gzipped else f[-1]
        file_type = False
        if suffix == 'bed' or suffix == 'vcf':
            file_type = suffix
        elif suffix == 'txt':
            if 'snps' in f:
                file_type = 'snpfile'
        if file_type:
            if verbose:
                print("File " + file + " type: " + file_type, file=sys.stderr)
            return(file_type, gzipped)
        else:
            if verbose:
                print("File " + file + " is not recognized", file=sys.stderr)
            return(False, False)

    #################################
    #     Create the dictionary     #
    #################################
    snps = {}

    if type(input) == str:
        input = [input]
    elif type(input) != list:
        raise TypeError('Argument to parse_snpfile must be a list')

    # Deal with directories
    if os.path.isdir(input[0]):
        if verbose:
            print('Directory detected', file=sys.stderr)
        if len(input) > 1:
            printerr(input[0] + ' is a directory, but you have provided too many arguments, exiting.')
            raise Exception('too many directories')
        files = [file for file in os.listdir(input[0]) if os.path.isfile(file)]
        for file in files:
            file_type, gzipped = name_parse(file)
            if file_type:
                snps.update(parse_file(file, file_type, gzipped))
    else:
        for file in input:
            file_type, gzipped = name_parse(file)
            if file_type:
                snps.update(parse_file(file, file_type, gzipped))

    return(snps)

def printerr(err):
    """ Print a message in red to STDERR """
    print('\033[91m' + err + '\x1b[0m', file=sys.stderr)

def log(output, logfile, printerr=False):
    """ Log to file with timestamp
        if printerr is True, will also print to stderr """
    import datetime

    timestamp   = datetime.datetime.now().strftime("%Y%m%d %H:%M:%S")
    output      = str(output)
    timeput     = ' | '.join([timestamp, output])

    print(timeput, file=logfile)
    if printerr:
        print(output, file=sys.stderr)

