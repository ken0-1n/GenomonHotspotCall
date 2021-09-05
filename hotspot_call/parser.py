#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: ken0-1n
"""

import sys
import argparse
from .version import __version__
from .run import hotspot_call_main

def create_parser():
    prog = "hotspotCall"
    parser = argparse.ArgumentParser(prog = prog)
    parser.add_argument("--version", action = "version", version = prog + "-" + __version__)
    subparsers = parser.add_subparsers()
                    
    def _create_mutation_parser(subparsers):
    
        mutation_parser = subparsers.add_parser("main", help = "calling hotspot")
        mutation_parser.add_argument("tumor_bam", metavar = "tumor.bam", type = str, help = "the path to the tumor bam file")
        mutation_parser.add_argument("control_bam", metavar = "control.bam", type = str, help = "the path to the control bam file")
        mutation_parser.add_argument("output_file", metavar = "output_file", type = str, help = "the path to the output file")
        mutation_parser.add_argument("hotspot_file", metavar = "hotspot_mutations.bed", type = str, help = "the bed format file that lists mutations")
        mutation_parser.add_argument('-O', '--print_format', choices = ['vcf','anno'], help = 'Print VCF or anno(TSV) format',  default = 'anno' )
        mutation_parser.add_argument('-S', metavar = 'samtools_params', type = str, default = "-B -q 20 -Q2 -d 10000000")
        mutation_parser.add_argument('-t', metavar = "min_tumor_misrate", default = 0.1, type = float)
        mutation_parser.add_argument('-c', metavar = "max_control_misrate", default = 0.1, type = float)
        mutation_parser.add_argument('-R', metavar = "T/N ratio_control", default = 0.1, type = float)
        mutation_parser.add_argument('-m', metavar = "min_lod_score", default = 8.0, type = float)
        mutation_parser.add_argument('-r', metavar = "rna.bam" , type = str, help = "the path to the RNA bam file")
        mutation_parser.add_argument( '-1', '--sample1', help = '1st sample name ( disease )', type = str, default = None)
        mutation_parser.add_argument( '-2', '--sample2', help = '2nd sample name ( control )', type = str, default = None)
        mutation_parser.add_argument( '-3', '--sample3', help = '3rd sample name ( rnaseq  )', type = str, default = None)
        mutation_parser.add_argument( '-f', '--ref_fa',  help = 'Reference genome', type = str, default = None)
        return mutation_parser
    
    mutation_parser = _create_mutation_parser(subparsers)
    mutation_parser.set_defaults(func = hotspot_call_main)
    return parser
    
