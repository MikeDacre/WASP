#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8 tabstop=4 expandtab shiftwidth=4 softtabstop=4
"""
================================================================================
                            WASP - file_handing

Simple functions for handing files, IO, and logging

 Last modified: 2015-03-12 16:11

================================================================================
"""
import pysam, sys

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

