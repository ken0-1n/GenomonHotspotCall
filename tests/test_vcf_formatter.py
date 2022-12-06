
import sys
import unittest, subprocess
import os, tempfile, shutil, filecmp
from hotspot_call import vcf_formatter as vf
from hotspot_call.fisher_info import FisherInfo


class TestVCFFormatter(unittest.TestCase):

    # get_ebpval():
    def test01(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 60
        fi.tumor_bases["a"] = 40
        fi.tumor_bases["C"] = 60
        fi.tumor_bases["c"] = 40
        fi.ctrl_bases["C"] = 30
        fi.ctrl_bases["c"] = 10
        fi.ctrl_bases["a"] = 1
        fi.chr = "chr1"
        fi.start = "10000"
        fi.end = "10001"
        fi.ref = "C"
        alt = "A"
        is_ctrl = True
        is_rna = False
        ret = vf.make_record(fi, alt, is_rna, is_ctrl)
        ans = "chr1\t10000\t.\tC\tA\t.\t.\tFP=9.021;LS=6.680\tDP:DPF:DPR:AD:ADF:ADR:AF:SB\t200:120:80:100:60:40:0.500:0.600\t41:30:11:1:0:1:0.024:0.000"
        self.assertTrue(ret == ans)
   
    def test02(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 60
        fi.tumor_bases["a"] = 40
        fi.tumor_bases["C"] = 60
        fi.tumor_bases["c"] = 40
        fi.ctrl_bases["C"] = 30
        fi.ctrl_bases["c"] = 10
        fi.ctrl_bases["a"] = 1
        fi.rna_bases["C"] = 20
        fi.rna_bases["c"] = 6
        fi.rna_bases["A"] = 10
        fi.rna_bases["a"] = 3
        fi.chr = "chr1"
        fi.start = "10000"
        fi.end = "10001"
        fi.ref = "C"
        alt = "A"
        is_ctrl = True
        is_rna = True
        ret = vf.make_record(fi, alt, is_rna, is_ctrl)
        ans = "chr1\t10000\t.\tC\tA\t.\t.\tFP=9.021;LS=6.680\tDP:DPF:DPR:AD:ADF:ADR:AF:SB\t200:120:80:100:60:40:0.500:0.600\t41:30:11:1:0:1:0.024:0.000\t39:30:9:13:10:3:0.333:0.769"
        self.assertTrue(ret == ans)
   
    def test03(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 60
        fi.tumor_bases["a"] = 40
        fi.tumor_bases["C"] = 60
        fi.tumor_bases["c"] = 40
        fi.ctrl_bases["C"] = 30
        fi.ctrl_bases["c"] = 10
        fi.ctrl_bases["a"] = 1
        fi.chr = "chr1"
        fi.start = "10000"
        fi.end = "10001"
        fi.ref = "C"
        alt = "A"
        is_ctrl = False
        is_rna = False
        ret = vf.make_record(fi, alt, is_rna, is_ctrl)
        ans = "chr1\t10000\t.\tC\tA\t.\t.\tB10=0.455;BM=0.500;B90=0.545;LS=6.680\tDP:DPF:DPR:AD:ADF:ADR:AF:SB\t200:120:80:100:60:40:0.500:0.600"
        self.assertTrue(ret == ans)
   
