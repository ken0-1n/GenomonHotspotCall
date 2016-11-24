
#! /usr/bin/env python
import math
import numpy
from scipy.stats import fisher_exact as fisher
import re

target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
remove_chr = re.compile( '\^.' )

class FisherInfo:

    def __init__(self):
    
        self.chr = 0
        self.start = 0
        self.end = 0
        self.ref = ""
        self.tumor_bases = {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                   }
        self.ctrl_bases = {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                   }
        self.rna_bases = {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                   }
        self.tumor_quals = {
                    "A": [],
                    "C": [],
                    "G": [],
                    "T": [],
                    "a": [],
                    "c": [],
                    "g": [],
                    "t": [] 
                   }
        self.ctrl_quals = {
                    "A": [],
                    "C": [],
                    "G": [],
                    "T": [],
                    "a": [],
                    "c": [],
                    "g": [],
                    "t": [] 
                   }
        self.rna_quals = {
                    "A": [],
                    "C": [],
                    "G": [],
                    "T": [],
                    "a": [],
                    "c": [],
                    "g": [],
                    "t": [] 
                   }


    def bases_format_process(self, read_bases, qual_list):
    
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

    def set_mpileup_data(self, mp_list):

        # Prepare mpileup data
        self.chr =  mp_list[0]
        self.start = mp_list[1]
        self.end = mp_list[1]

        tumor_bases = self.bases_format_process(mp_list[4], mp_list[5])
        for base in tumor_bases:
            self.add_tumor_base(base) 

        for base, qual in zip(tumor_bases, mp_list[5]):
            self.add_tumor_quals(base, qual)

        ctrl_bases = self.bases_format_process(mp_list[7], mp_list[8])
        for base in ctrl_bases:
            self.add_ctrl_base(base)

        for base, qual in zip(ctrl_bases, mp_list[8]):
            self.add_ctrl_quals(base, qual)

        if len(mp_list) > 10:
            rna_bases = self.bases_format_process(mp_list[10], mp_list[11])
            for base in rna_bases:
                self.add_rna_base(base)

            for base, qual in zip(rna_bases, mp_list[11]):
                self.add_rna_quals(base, qual)


    def set_ref(self,ref):
        self.ref = ref

    def add_base(self,bases,base):
        if base in 'ATGCatgc':
            bases[base] += 1

    def add_tumor_base(self, base):
        self.add_base(self.tumor_bases, base)

    def add_ctrl_base(self, base):
        self.add_base(self.ctrl_bases, base)

    def add_rna_base(self, base):
        self.add_base(self.rna_bases, base)

    def add_quals(self, quals, base, qual):
        if base in 'ATGCatgc':
            ord_qual = (int(ord(qual))-33)
            q = quals[base]
            q.append(min(ord_qual,41))

    def add_tumor_quals(self, base, qual):
        self.add_quals(self.tumor_quals, base, qual)

    def add_ctrl_quals(self, base, qual):
        self.add_quals(self.ctrl_quals, base, qual)

    def add_rna_quals(self, base, qual):
        self.add_quals(self.rna_quals, base, qual)

    def get_depth(self, bases):
        count = 0
        for n in "ACGTacgt":
            count += bases[n]
        return count

    def get_tumor_depth(self):
        return self.get_depth(self.tumor_bases)

    def get_ctrl_depth(self):
        return self.get_depth(self.ctrl_bases)

    def get_rna_depth(self):
        return self.get_depth(self.rna_bases)

    def get_depth_plus_strand(self, bases):
        count = 0
        for n in "ACGT":
            count += bases[n]
        return count

    def get_tumor_depth_plus_strand(self):
        return self.get_depth_plus_strand(self.tumor_bases)

    def get_ctrl_depth_plus_strand(self):
        return self.get_depth_plus_strand(self.ctrl_bases)

    def get_rna_depth_plus_strand(self):
        return self.get_depth_plus_strand(self.rna_bases)

    def get_depth_minus_strand(self, bases):
        count = 0
        for n in "acgt":
            count += bases[n]
        return count

    def get_tumor_depth_minus_strand(self):
        return self.get_depth_minus_strand(self.tumor_bases)

    def get_ctrl_depth_minus_strand(self):
        return self.get_depth_minus_strand(self.ctrl_bases)

    def get_rna_depth_minus_strand(self):
        return self.get_depth_minus_strand(self.rna_bases)

    def get_tumor_base_total(self, base):
        return (self.tumor_bases[base.upper()] + self.tumor_bases[base.lower()])

    def get_ctrl_base_total(self, base):
        return (self.ctrl_bases[base.upper()] + self.ctrl_bases[base.lower()])

    def get_rna_base_total(self, base):
        return (self.rna_bases[base.upper()] + self.rna_bases[base.lower()])

    def get_tumor_base_plus_strand(self, base):
        return (self.tumor_bases[base.upper()])

    def get_ctrl_base_plus_strand(self, base):
        return (self.ctrl_bases[base.upper()])

    def get_rna_base_plus_strand(self, base):
        return (self.rna_bases[base.upper()])

    def get_tumor_base_minus_strand(self, base):
        return (self.tumor_bases[base.lower()])

    def get_ctrl_base_minus_strand(self, base):
        return (self.ctrl_bases[base.lower()])

    def get_rna_base_minus_strand(self, base):
        return (self.rna_bases[base.lower()])

    def get_misrate(self,mis_base_count,depth):
        if mis_base_count == 0:
            return float(0)
        else:
           return (mis_base_count / float(depth))

    def get_tumor_misrate(self,base):
        return self.get_misrate(self.get_tumor_base_total(base), self.get_tumor_depth())

    def get_ctrl_misrate(self,base):
        return self.get_misrate(self.get_ctrl_base_total(base), self.get_ctrl_depth())

    def get_rna_misrate(self,base):
        return self.get_misrate(self.get_rna_base_total(base), self.get_rna_depth())

    def get_strand_ratio(self,mis_base_count_plus,mis_base_count_minus):
        if (mis_base_count_plus + mis_base_count_minus) == 0:
            return float(-1)
        elif mis_base_count_plus == 0:
            return float(0)
        else:
            return (mis_base_count_plus / float(mis_base_count_plus + mis_base_count_minus))
            
    def get_tumor_strand_ratio(self,base):
        return self.get_strand_ratio(self.get_tumor_plus_strand(base), self.get_tumor_minus_strand(base))
            
    def get_ctrl_strand_ratio(self, base):
        return self.get_strand_ratio(self.get_ctrl_plus_strand(base), self.get_ctrl_minus_strand(base))

    def get_rna_strand_ratio(self, base):
        return self.get_strand_ratio(self.get_rna_plus_strand(base), self.get_rna_minus_strand(base))

    def get_fisher_pvalue(self,base):
        odds_ratio, fisher_pvalue = fisher(
        ((int(self.get_tumor_base_total(self.ref)), int(self.get_ctrl_base_total(self.ref))),
         (int(self.get_tumor_base_total(base)),  int(self.get_ctrl_base_total(base)))),
          alternative='two-sided'
        )
        val = float(0.0)
        if fisher_pvalue < 10**(-60):
            val = float(60.0)
        elif fisher_pvalue  > 1.0 - 10**(-10) :
            val = float(0.0)
        else:
            val = -math.log( fisher_pvalue, 10 )
        return val

    def lod_qual(self, base):
        score = float(0)
        for qual in self.tumor_quals[base]:
            q = float(qual)
            p = 10**-(q/10)
            score += -math.log(p/(1-p),10)
        return score

    def get_lod_score(self,base):
        return (self.lod_qual(base.upper()) + self.lod_qual(base.lower()))

    def get_lod_score_plus_strand(self,base):
        return self.lod_qual(base.upper())
    
    def get_lod_score_minus_strand(self,base):
        return self.lod_qual(base.lower())

    def get_score_median(self,base):
        med = 0
        if len(self.tumor_quals[base]) != 0 or len(self.tumor_quals[base.lower()]) != 0: 
            alt_array = self.tumor_quals[base] + self.tumor_quals[base.lower()]
            med = numpy.median(alt_array)
        return med

