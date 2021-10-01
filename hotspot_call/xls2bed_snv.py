import gzip, os, sys
import pandas as pd
import subprocess

def make_hotspot_dict(in_tmp_hotspot):

    d_symbol_aa = {}    # for SNV
    d_symbol_codon = {} # for MNV
    d_symbol_pos = {}   # for SNV and MNV

    with open(in_tmp_hotspot,"r") as hin:
        for line in hin:
            if line.startswith("Hugo_Symbol"): continue
            F = line.split("\t")

            symbol,posaa=F[0],F[1]
            refaa = F[4].split(":")[0].replace("splice","")
            altaa = F[8].split(":")[0].replace("splice","")

            # if silent mutation
            if refaa == altaa and posaa.find('splice') == -1: altaa = "="
            # if codon pos is 1
            if posaa == "1" and posaa.find('splice') == -1: altaa = "?"

            aa_change = "p."+refaa+posaa+altaa
            d_symbol_aa.setdefault(symbol, set()).add(aa_change)
    
            for genome_pos in F[10].split("|"):
                d_symbol_pos.setdefault(symbol, set()).add(genome_pos.split("_")[0])
    
            for codon in F[9].split("|"):
                d_symbol_codon.setdefault(symbol, set()).add(codon.split(":")[0])

    return d_symbol_aa, d_symbol_codon, d_symbol_pos
    

def get_hotspot_record(in_maf, d_symbol_aa, d_symbol_codon, d_symbol_pos):

    d_mut = {}
    with gzip.open(in_maf,"rt") as hin:
        for line in hin:
            if line.startswith("#"): continue
            if line.startswith("Hugo_Symbol"): continue

            F = line.rstrip("\n").split("\t")
            if F[9] not in ["SNP","DNP","ONP","TNP"]: continue

            sample = F[15]
            if sample.startswith("P-"):
                sample = "P-"+sample.split("-")[1]
            hgvsc, hgvsp_short = F[34], F[36]
            codons = F[55]
   
            # error check 1
            if len(set([F[10],F[11],F[12]])) > 2 and F[11]!= "NA":
                print(F[10]+","+F[11]+","+F[12], file=sys.stderr)
                sys.exit(1)
    
            # error check 2
            if F[10] == F[12]:
                print(F[10]+","+F[11]+","+F[12], file=sys.stderr)
                sys.exit(1)
    
            flag = False
            # SNP
            if (F[0] in d_symbol_aa and hgvsp_short in d_symbol_aa[F[0]]):
                if F[4]+":"+F[5] in d_symbol_pos[F[0]]:
                    flag = True
            # MNP
            elif F[9] != "SNP" and F[0] in d_symbol_codon and codons in d_symbol_codon[F[0]]:
                if F[4]+":"+F[5] in d_symbol_pos[F[0]]:
                    flag = True
    
            if flag:
                # key = symbol,chrom,start,end,ref,alt2,sample
                sample_key = "\t".join([F[0],F[4],F[5],F[6],F[10],F[12],sample])
                if sample_key in d_mut: continue
                d_mut[sample_key] = "\t".join([F[4],F[5],F[6],F[10],F[12],F[0],hgvsc,hgvsp_short,F[9]])

    return d_mut


def database_snv_main(args):

    out_prefix, ext = os.path.splitext(args.out_snv_database)
    in_tmp_hotspot = out_prefix +".tsv"

    # xls to tsv
    input_book = pd.ExcelFile(args.in_xls_hotspot)
    input_sheet_name = input_book.sheet_names
    input_sheet_df = input_book.parse(input_sheet_name[0])
    input_sheet_df.to_csv(in_tmp_hotspot, sep='\t', index=False)

    d_hotspot, d_codon, d_genome_pos_symbol = make_hotspot_dict(in_tmp_hotspot)
    d_mut = get_hotspot_record(args.in_maf_hotspot, d_hotspot, d_codon, d_genome_pos_symbol)
   
    with open(out_prefix +".bed", 'w') as hout: 
        for key in d_mut:
           l = d_mut[key].split('\t')
           print("\t".join(["chr"+l[0],str(int(l[1])-1),l[2],l[3],l[4],l[5],l[6],l[7],l[8]]), file=hout)

    subprocess.check_call(['sort', '-k', '1,1', '-k', '2,2n', '-k', '3,3n', '-k', '4,4', '-k', '5,5', '-o', out_prefix+".sorted.bed", out_prefix+".bed"]) 

    subprocess.check_call(['uniq', out_prefix+".sorted.bed", out_prefix+".uniq.bed"]) 

    subprocess.check_call(['liftOver', '-bedPlus=3', out_prefix+".uniq.bed", args.map_chain, out_prefix+".liftedover.bed", out_prefix+".unmapped.bed"]) 

    with open(out_prefix+".liftedover.bed", 'r') as hin, open(args.out_snv_database, 'w') as hout: 
        for line in hin:
           l = line.rstrip('\n').split('\t')
           print("\t".join([l[0],str(int(l[1])+1),l[2],l[3],l[4],l[5],l[6],l[7],l[8]]), file=hout)

    if not args.debug:
        os.remove(out_prefix +".bed")
        os.remove(out_prefix +".tsv")
        os.remove(out_prefix+".sorted.bed")
        os.remove(out_prefix+".uniq.bed")
        os.remove(out_prefix+".liftedover.bed")

