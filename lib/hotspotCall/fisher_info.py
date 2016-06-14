
#! /usr/bin/env python
import math
from scipy.stats import fisher_exact as fisher

class FisherInfo:

    def __init__(self):
    
        self.chr = 0
        self.start = 0
        self.end = 0
        self.ref = ""
        self.target_base_num = {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                   }
        self.control_base_num = {
                    "A": 0,
                    "C": 0,
                    "G": 0,
                    "T": 0,
                    "a": 0,
                    "c": 0,
                    "g": 0,
                    "t": 0
                   }

    def add_target_base(self, base):
        if base in 'ATGCatgc':
            self.target_base_num[base] += 1

    def add_control_base(self, base):
        if base in 'ATGCatgc':
            self.control_base_num[base] += 1
    
    def get_target_depth_total(self):
        count = 0
        for n in "ACGTacgt":
            count += self.target_base_num[n]
        return count

    def get_control_depth_total(self):
        count = 0
        for n in "ACGTacgt":
            count += self.control_base_num[n]
        return count

    def get_target_depth_plus_strand(self):
        count = 0
        for n in "ACGT":
            count += self.target_base_num[n]
        return count

    def get_control_depth_plus_strand(self):
        count = 0
        for n in "ACGT":
            count += self.control_base_num[n]
        return count

    def get_target_depth_minus_strand(self):
        count = 0
        for n in "acgt":
            count += self.target_base_num[n]
        return count

    def get_control_depth_minus_strand(self):
        count = 0
        for n in "acgt":
            count += self.control_base_num[n]
        return count

    def get_target_total(self, base):
        return (self.target_base_num[base.upper()] + self.target_base_num[base.lower()])

    def get_control_total(self, base):
        return (self.control_base_num[base.upper()] + self.control_base_num[base.lower()])

    def get_target_plus_strand(self, base):
        return (self.target_base_num[base.upper()])

    def get_control_plus_strand(self, base):
        return (self.control_base_num[base.upper()])

    def get_target_minus_strand(self, base):
        return (self.target_base_num[base.lower()])

    def get_control_minus_strand(self, base):
        return (self.control_base_num[base.lower()])

    def get_target_A_total(self):
        return (self.target_base_num["A".upper()] + self.target_base_num["A".lower()])

    def get_target_C_total(self):
        return (self.target_base_num["C".upper()] + self.target_base_num["C".lower()])

    def get_target_G_total(self):
        return (self.target_base_num["G".upper()] + self.target_base_num["G".lower()])

    def get_target_T_total(self):
        return (self.target_base_num["T".upper()] + self.target_base_num["T".lower()])

    def get_control_A_total(self):
        return (self.control_base_num["A".upper()] + self.control_base_num["A".lower()])

    def get_control_C_total(self):
        return (self.control_base_num["C".upper()] + self.control_base_num["C".lower()])

    def get_control_G_total(self):
        return (self.control_base_num["G".upper()] + self.control_base_num["G".lower()])

    def get_control_T_total(self):
        return (self.control_base_num["T".upper()] + self.control_base_num["T".lower()])

    def get_target_misrate(self,base):
        mis_base_count = self.get_target_total(base)
        if mis_base_count == 0:
            return float(0)
        else:
           return (mis_base_count / float(self.get_target_depth_total()))

    def get_control_misrate(self,base):
        mis_base_count = self.get_control_total(base)
        if mis_base_count == 0:
            return float(0)
        else:
           return (mis_base_count / float(self.get_control_depth_total()))

    def get_target_strand_ratio(self,base):
        mis_base_count_plus = self.get_target_plus_strand(base)
        mis_base_count_minus = self.get_target_minus_strand(base)
        if (mis_base_count_plus + mis_base_count_minus) == 0:
            return float(-1)
        elif mis_base_count_plus == 0:
            return float(0)
        else:
            return (mis_base_count_plus / float(mis_base_count_plus + mis_base_count_minus))
            
    def get_control_strand_ratio(self, base):
        mis_base_count_plus = self.get_control_plus_strand(base)
        mis_base_count_minus = self.get_control_minus_strand(base)
        if (mis_base_count_plus + mis_base_count_minus) == 0:
            return float(-1)
        elif mis_base_count_plus == 0:
            return float(0)
        else:
            return (mis_base_count_plus / float(mis_base_count_plus + mis_base_count_minus))

    def get_fisher_pvalue(self,base):
        odds_ratio, fisher_pvalue = fisher(
        ((int(self.get_target_total(self.ref)),
          int(self.get_control_total(self.ref))),
         (int(self.get_target_total(base)),
          int(self.get_control_total(base)))),
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


