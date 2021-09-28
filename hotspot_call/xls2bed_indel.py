import gzip, os, sys
import pandas as pd
import re

def make_hotspot_dict(in_tmp_hotspot):

    d_hotspot = {}
    d_genome_pos = {}
    with open(in_tmp_hotspot,"r") as hin:
        for line in hin:
            if line.startswith("Hugo_Symbol"): continue
            F = line.split("\t")
            symbol,posaa=F[0],F[1]
            aa_change = "p."+F[9].split(":")[0]
            d_hotspot.setdefault(symbol, []).append(aa_change)
    
            for genome_pos in F[11].split("|"):
                genome_pos = genome_pos.split("_")[0]
                d_genome_pos.setdefault(symbol+"\t"+aa_change, []).append(genome_pos)

    return d_hotspot, d_genome_pos 


def get_hotspot_record(in_hotspot, d_hotspot, d_genome_pos):
    l_mut = []
    with gzip.open(in_hotspot,"rt") as hin:
        for line in hin:
            if line.startswith("#"): continue
            if line.startswith("Hugo_Symbol"): continue
            line = line.rstrip("\n")
            F = line.split("\t")
            symbol = F[0]
            chrom, start, end = F[4],F[5],F[6] 
            mut_type = F[9]
            ref, alt1, alt2 = F[10],F[11],F[12] 
            hgvsc, hgvsp_short = F[34], F[36]
   
            # check1
            if len(set([ref,alt1,alt2])) > 2 and alt1 != "NA":
                print(ref+","+alt1+","+alt2, file=sys.stderr)
                sys.exit(1)
    
            # check2
            if alt2 == ref:
                print(ref+","+alt1+","+alt2, file=sys.stderr)
                sys.exit(1)
    
            if symbol in d_hotspot and hgvsp_short in d_hotspot[symbol]:
                if chrom+":"+start in d_genome_pos[symbol+"\t"+hgvsp_short]:
                    l_mut.append("\t".join([chrom,start,end,ref,alt2,symbol,hgvsc,hgvsp_short,mut_type]))

    return l_mut


def summary_hotspot_record(out_indel, l_mut):

    with open(out_indel, 'w') as hout_indel: 
    
        for val in list(set(l_mut)):
           mut_type = val.split("\t")[8]
           if mut_type == "INS" or mut_type == "DEL": 
              print(val, file=hout_indel)
        

def database_indel_main(args):

    in_xls_hotspot_name, ext = os.path.splitext(args.in_xls_hotspot)
    in_tmp_hotspot = in_xls_hotspot_name +".tsv"

    # xls to tsv
    input_book = pd.ExcelFile(args.in_xls_hotspot)
    input_sheet_name = input_book.sheet_names
    input_sheet_df = input_book.parse(input_sheet_name[1])
    input_sheet_df.to_csv(in_tmp_hotspot, sep='\t', index=False)

    d_hotspot, d_genome_pos = make_hotspot_dict(in_tmp_hotspot)
    l_mut = get_hotspot_record(args.in_maf_hotspot, d_hotspot, d_genome_pos)
    summary_hotspot_record(args.out_indel_database, l_mut)

    os.remove(in_tmp_hotspot)


