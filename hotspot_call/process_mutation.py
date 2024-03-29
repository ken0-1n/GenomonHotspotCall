#! /usr/bin/env python
import sys, os, re, subprocess
from . import anno_formatter
from . import vcf_formatter
from .fisher_info import FisherInfo

def read_hotspot_file(hotspot_file):
    tmp_list = []
    hIN = open(hotspot_file, 'r')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if len(F[3]) == 1 and len(F[4]) == 1 and F[3] in 'ACGT' and F[4] in 'ACGT':
            key = F[0] +"\t"+ F[1] +"\t"+ F[2] +"\t"+ F[3] +"\t"+ F[4]
            tmp_list.append(key)
    hIN.close()

    tmp_sorted_list = sorted(set(tmp_list), key=tmp_list.index)
    
    tmp_dict = {}
    for line in tmp_sorted_list:
        F = line.split('\t')
        key = F[0] +"\t"+ F[1] +"\t"+ F[2] +"\t"+ F[3]
        if key in tmp_dict:
            var = tmp_dict[key]
            tmp_dict[key] = var +","+ F[4]
        else:
            tmp_dict[key] = F[4]

    return tmp_dict


def print_anno_header(is_rna, hOUT, is_ctrl):

    header_str = "Chr\tStart\tEnd\tRef\tAlt"
    if is_ctrl:
        header_str = header_str +"\tdepth_tumor\tvariantNum_tumor\tdepth_normal\tvariantNum_normal\tbases_tumor\tbases_normal\tA_C_G_T_tumor\tA_C_G_T_normal\tmisRate_tumor\tstrandRatio_tumor\tmisRate_normal\tstrandRatio_normal\tP-value(fisher)\tscore"
    else:
        header_str = header_str +"\tdepth_tumor\tvariantNum_tumor\tbases_tumor\tbases_normal\tA_C_G_T_tumor\tmisRate_tumor\tstrandRatio_tumor\t10%_posterior_quantile\tposterior_mean\t90%_posterior_quantile\tscore"

    if is_rna:
        header_str = header_str +"\tdepth_RNA\tvariant_RNA\tbases_RNA\ttmisRate_RNA"
    print(header_str, file=hOUT)


def print_vcf_header(is_rna, sample1, sample2, sample_rna, ref_fa, ref_dict, hOUT, is_ctrl):

    print('##fileformat=VCFv4.2',file=hOUT)
    # print info and format
    if is_ctrl:
        print('##INFO=<ID=FP,Number=1,Type=Float,Description="Minus logarithm of the p-value by Fishers exact test">',file=hOUT)
    else:
        print('##INFO=<ID=B10,Number=1,Type=Float,Description="10% posterior quantile of the beta distribution">',file=hOUT)
        print('##INFO=<ID=BM,Number=1,Type=Float,Description="Posterior mean">',file=hOUT)
        print('##INFO=<ID=B90,Number=1,Type=Float,Description="90% posterior quantile of the beta distribution">',file=hOUT)
    print('##INFO=<ID=LS,Number=1,Type=Float,Description="LOD Score of Hotspot Call">',file=hOUT)
    print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',file=hOUT)
    print('##FORMAT=<ID=DPF,Number=1,Type=Integer,Description="Read depth in the forward strand">',file=hOUT)
    print('##FORMAT=<ID=DPR,Number=1,Type=Integer,Description="Read depth in the reverse strand">',file=hOUT)
    print('##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allelic depth">',file=hOUT)
    print('##FORMAT=<ID=ADF,Number=1,Type=Integer,Description="Allelic depth in the forward strand">',file=hOUT)
    print('##FORMAT=<ID=ADR,Number=1,Type=Integer,Description="Allelic depth in the reverse strand">',file=hOUT)
    print('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">',file=hOUT)
    print('##FORMAT=<ID=SB,Number=1,Type=Float,Description="Strand bias">',file=hOUT)

    # print reference information
    hIN = open(ref_dict)
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if F[0] == '@SQ':
            ID = F[1].replace('SN:','')
            length = F[2].replace('LN:','')
            print('##contig=<ID='+ID+',length='+length+'>',file=hOUT)
    print('##reference='+ref_fa,file=hOUT)
    hIN.close()

    # print_header
    samples = sample1
    if is_ctrl:
        samples = samples +"\t"+ sample2
    if is_rna:
        samples = samples +"\t"+ sample_rna
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+samples,file=hOUT)


def make_pileup_cmmand(mutReg, m_params, is_rna, is_ctrl, bam_tumor, bam_control, bam_rna):
    mpileup_cmd = ["samtools", "mpileup", "-r", mutReg]
    mpileup_cmd.extend(m_params)

    if is_rna:
        if is_ctrl:
            mpileup_cmd.extend([bam_tumor, bam_control, bam_rna])
        else:
            print('Error: This tool not supported without entering the matced control and rnaseq bam.', file=sys.stderr)
            sys.exit(1)
    else:
        if is_ctrl:
            mpileup_cmd.extend([bam_tumor, bam_control])
        else:
            mpileup_cmd.extend([bam_tumor])
    return mpileup_cmd
        

def is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl):
    ret = False
    if fi.get_lod_score(alt) < min_lod_score or fi.get_tumor_misrate(alt) < min_tumor_misrate:
         ret = True
    if is_ctrl:
        if fi.get_ctrl_misrate(alt) > max_ctrl_misrate or fi.get_ctrl_misrate(alt) > (fi.get_tumor_misrate(alt) * ratio_ctrl):
             ret = True
    return ret


def call(hotspot_file, output_file, bam_tumor, bam_control, mpileup_params, min_tumor_misrate, max_ctrl_misrate, bam_rna, min_lod_score, ratio_ctrl, is_anno, sample1, sample2, sample_rna, ref_fa):

    is_rna = True if bam_rna else False
    is_ctrl = True if bam_control else False
    hOUT = open(output_file, 'w')
    FNULL = open(os.devnull, 'w')

    if is_anno:
        print_anno_header(is_rna, hOUT, is_ctrl)
    else:
        ref_name, ext = os.path.splitext(ref_fa)
        print_vcf_header(is_rna, sample1, sample2, sample_rna, ref_fa, ref_name+".dict", hOUT, is_ctrl)

    m_params = mpileup_params.split(" ")
    hotspot_dict = read_hotspot_file(hotspot_file)
    for key in hotspot_dict:

        F = key.split('\t')
        mutReg = F[0] +":"+ F[1] +"-"+ F[2]
        hotspot_ref = F[3]
        hotspot_alts = hotspot_dict[key].split(',')

        # TODO: error message
        if F[1] != F[2]:
            print("Invalid position in the hotspot database: "+ F[0] +"\t"+ F[1] +"\t"+ F[2], file=sys.stderr)
            continue
        # TODO: error message
        for tmp_alt in hotspot_alts:
            if tmp_alt not in "ACGTacgt":
                print("Invalid Alt in the mutations.bed: "+ F[0] +"\t"+ F[1] +"\t"+ F[2] +"\t"+ hotspot_alts, file=sys.stderr)
           
        mpileup_cmd = make_pileup_cmmand(mutReg, m_params, is_rna, is_ctrl, bam_tumor, bam_control, bam_rna)

        # print mpileup_cmd
        pileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr = FNULL)
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            fi = FisherInfo()
            mp_list = mpileup.decode().strip('\n').split( '\t' )
            fi.set_ref(hotspot_ref)
            fi.set_mpileup_data(mp_list)
            for alt in hotspot_alts:
                if is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl): continue
                if is_anno:
                    record = anno_formatter.make_record(fi, alt, is_rna, is_ctrl)
                else:
                    record = vcf_formatter.make_record(fi, alt, is_rna, is_ctrl)
                print(record,file=hOUT)

    FNULL.close()
    hOUT.close()

