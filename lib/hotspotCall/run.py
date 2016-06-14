#! /usr/bin/env python

# import process_vcf
import process_mutation
# import get_eb_score
# import sys, os, subprocess, math, re, multiprocessing 
import sys, os, subprocess
# import vcf, pysam, numpy
# import vcf


def worker_vcf(targetBamPath, controlBamPath, outputPath, hotspotList, mapping_qual_thres, base_qual_thres, min_allele_thres):

    # generate file
    process_mutation.call(hotspotList, outputPath + '.target.vcf', targetBamPath, controlBamPath, mapping_qual_thres, base_qual_thres, min_allele_thres)


def worker_anno(targetBamPath, controlBamPath, outputPath, hotspotList, mapping_qual_thres, base_qual_thres, min_allele_thres):

    # generate file
    process_mutation.call(hotspotList, outputPath + '.target.anno', targetBamPath, controlBamPath, mapping_qual_thres, base_qual_thres, min_allele_thres)


def main(args):

    # should add validity check for arguments
    targetBamPath = args.targetBamPath
    controlBamPath = args.controlBamPath
    outputPath = args.outputPath
    hotspotList= args.hotspotList

    mapping_qual_thres = args.q
    base_qual_thres = args.Q
    min_allele_thres = args.m
    is_anno = True if args.f == 'anno' else False

    # file existence check
    if not os.path.exists(targetBamPath):
        print >> sys.stderr, "No target bam file: " + targetBamPath
        sys.exit(1)

    if not os.path.exists(targetBamPath + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", targetBamPath)):
        print >> sys.stderr, "No index for target bam file: " + targetBamPath
        sys.exit(1)

    if not os.path.exists(controlBamPath):
        print >> sys.stderr, "No control bam file: " + controlBamPath
        sys.exit(1)

    if not os.path.exists(controlBamPath + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", controlBamPath)):
        print >> sys.stderr, "No index for control bam file: " + controlBamPath
        sys.exit(1)

    if not os.path.exists(hotspotList):
        print >> sys.stderr, "No hotspot mutations list: " + hotspotList
        sys.exit(1)

    outputDir = os.path.dirname(outputPath)
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    # non multi-threading mode
    if is_anno == True:
        worker_anno(targetBamPath, controlBamPath, outputPath, hotspotList, mapping_qual_thres, base_qual_thres, min_allele_thres)
    else: 
        worker_vcf(targetBamPath, controlBamPath, outputPath, hotspotList, mapping_qual_thres, base_qual_thres, min_allele_thres)

    # delete intermediate files
    # os.unlink(outputPath + '.target.pileup')


