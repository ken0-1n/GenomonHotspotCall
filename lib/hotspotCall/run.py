#! /usr/bin/env python

import process_mutation
import sys, os, subprocess

def main(args):

    # should add validity check for arguments
    tumor_bam = args.tumor_bam
    control_bam = args.control_bam
    output = args.output_file
    hotspot_file= args.hotspot_file
    # options
    mpileup_params = args.S
    min_tumor_misrate = args.t
    max_ctrl_misrate = args.c
    ratio_control = args.R
    min_lod_score = args.m
    rna_bam = args.r
    is_rna = True if rna_bam else False

    # file existence check
    if not os.path.exists(tumor_bam):
        print >> sys.stderr, "No tumor bam file: " + tumor_bam
        sys.exit(1)

    if not os.path.exists(tumor_bam + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", tumor_bam)):
        print >> sys.stderr, "No index for tumor bam file: " + tumor_bam
        sys.exit(1)

    if not os.path.exists(control_bam):
        print >> sys.stderr, "No control bam file: " + control_bam
        sys.exit(1)

    if not os.path.exists(control_bam + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", control_bam)):
        print >> sys.stderr, "No index for control bam file: " + control_bam
        sys.exit(1)

    if not os.path.exists(hotspot_file):
        print >> sys.stderr, "No hotspot mutations list: " + hotspot_file
        sys.exit(1)

    if is_rna:
        if not os.path.exists(rna_bam):
            print >> sys.stderr, "No rna bam file: " + rna_bam
            sys.exit(1)

        if not os.path.exists(rna_bam + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", rna_bam)):
            print >> sys.stderr, "No index for rna bam file: " + rna_bam
            sys.exit(1)

    process_mutation.call(hotspot_file, output, tumor_bam, control_bam, mpileup_params, min_tumor_misrate, max_ctrl_misrate, rna_bam, min_lod_score, ratio_control)

