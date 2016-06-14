#! /usr/bin/env python
import sys, os, re, subprocess
import anno_formatter
from fisher_info import FisherInfo

target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
remove_chr = re.compile( '\^.' )

def bases_format_process(read_bases, qual_list):

    deleted = 0
    iter = target.finditer( read_bases )
    for m in iter:
        site = m.start()
        type = m.group( 1 )
        num = m.group( 2 )
        bases = m.group( 3 )[ 0:int( num ) ]
        read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int( num ) + len( num ) + 1 - deleted: ]
        deleted += 1 + len( num ) + int( num )

    # Remove '^.' and '$'
    read_bases = remove_chr.sub( '', read_bases )
    read_bases = read_bases.translate( None, '$' ) 

    # Error check
    if len( read_bases ) != len( qual_list ):
        print >> sys.stderr, ("mpileup data is not good: {0}, {1}".format( read_bases, read_bases ))
        return None

    # Count mismatch
    return read_bases

def count_pileup(mpileup, hotspot_ref, hotspot_alts):

    # Prepare mpileup data
    mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
    # init data instance 
    fi = FisherInfo()
    fi.chr =  mp_list[0]
    fi.start = mp_list[1]
    fi.end = mp_list[1]
    fi.ref = hotspot_ref.upper()
    alts = []
    for alt in hotspot_alts:
        if not alt.upper() in alts:
            alts.append(alt.upper())
    fi.alts = alts

    target_bases = bases_format_process(mp_list[4], mp_list[5])
    for base in target_bases:
         fi.add_target_base(base) 
    
    control_bases = bases_format_process(mp_list[7], mp_list[8])
    for base in control_bases:
         fi.add_control_base(base)

    return fi

def read_mutation_file(mutation_list):
    tmp_list = []
    hIN = open(mutation_list, 'r')
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


def call(mutation_list, output_file, bam_tumor, bam_control, mq_thres, bq_thres, min_allele_freq):

    hOUT = open(output_file, 'w')
    FNULL = open(os.devnull, 'w')

    header_str = "#chr\tstart\tend\tref\talt\tdepth_tumor\tvariantNum_tumor\tdepth_normal\tvariantNum_normal\tbases_tumor\tbases_normal\tA,C,G,T_tumor\tA,C,G,T_normal\tmisRate_tumor\tstrandRatio_tumor\tmisRate_normal\tstrandRatio_normal\tP-value(fisher)"
    print >> hOUT, header_str

    mutation_dict = read_mutation_file(mutation_list)
    for key in mutation_dict:

        F = key.split('\t')
        mutReg = F[0] +":"+ F[1] +"-"+ F[2]
        hotspot_ref = F[3]
        hotspot_alts = mutation_dict[key].split(',')

        # TODO: error message
        if F[1] != F[2]: continue
        # TODO: error message
        for tmp_alt in hotspot_alts:
            if tmp_alt not in "ACGTacgt":
                print >> sys.stderr, "Invalid Alt in the mutations.bed: "+ F[0] +"\t"+ F[1] +"\t"+ hotspot_alts
              
        mpileup_cmd = ["samtools", "mpileup", "-B", "-d", "10000000", "-q", str(mq_thres), "-Q", str(bq_thres), "-r", mutReg, bam_tumor, bam_control]
        pileup = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE, stderr = FNULL)
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            fi = count_pileup(mpileup, hotspot_ref, hotspot_alts)
            for alt in fi.alts:
                if fi.get_target_misrate(alt) < min_allele_freq: continue
                record = anno_formatter.make_record(fi, alt)
                print >> hOUT, record

    FNULL.close()
    hOUT.close()

