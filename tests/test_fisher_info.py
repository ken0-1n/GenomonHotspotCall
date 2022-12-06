
import sys
import unittest, subprocess
import os, tempfile, shutil, filecmp
from hotspot_call.fisher_info import FisherInfo


class TestFisherInfo(unittest.TestCase):

    def test01(self):
        fi = FisherInfo()
        read_bases = "$GGGGGAAAAAggaa"
        qual_list = "IIIIIIIIIIIIII"
        ret = fi.bases_format_process(read_bases, qual_list)
        ans = "GGGGGAAAAAggaa"
        self.assertTrue(ret == ans)

    def test02(self):
        fi = FisherInfo()
        fi.ref = "G"
        read_bases = "$.....AAAAA,,aa"
        qual_list = "IIIIIIIIIIIIII"
        ret = fi.bases_format_process(read_bases, qual_list)
        ans = "GGGGGAAAAAggaa"
        self.assertTrue(ret == ans)
        
    def test03(self):
        fi = FisherInfo()
        fi.ref = "G"
        read_bases = "$.....-5CCCCCAAAAA,,aa"
        qual_list = "IIIIIIIIIIIIII"
        ret = fi.bases_format_process(read_bases, qual_list)
        ans = "GGGGGAAAAAggaa"
        self.assertTrue(ret == ans)
        
    def test11(self):
        fi = FisherInfo()
        fi.ref = "G"
        mp_list = []
        mp_list.append("chr1")
        mp_list.append("10000")
        mp_list.append("10001")
        mp_list.append("")
        mp_list.append("$.....AAAA,,a")
        mp_list.append("IIIIIIIIIIII")
        fi.set_mpileup_data(mp_list)
        self.assertTrue(fi.tumor_bases["G"] == 5)
        self.assertTrue(fi.tumor_bases["g"] == 2)
        self.assertTrue(fi.tumor_bases["A"] == 4)
        self.assertTrue(fi.tumor_bases["a"] == 1)
        
    def test12(self):
        fi = FisherInfo()
        fi.ref = "G"
        mp_list = []
        mp_list.append("chr1")
        mp_list.append("10000")
        mp_list.append("10001")
        mp_list.append("")
        mp_list.append("$.....AAAA,,a")
        mp_list.append("IIIIIIIIIIII")
        mp_list.append("")
        mp_list.append("......AA")
        mp_list.append("IIIIIIII")
        fi.set_mpileup_data(mp_list)
        self.assertTrue(fi.tumor_bases["G"] == 5)
        self.assertTrue(fi.tumor_bases["g"] == 2)
        self.assertTrue(fi.tumor_bases["A"] == 4)
        self.assertTrue(fi.tumor_bases["a"] == 1)
        self.assertTrue(fi.ctrl_bases["G"] == 6)
        self.assertTrue(fi.ctrl_bases["g"] == 0)
        self.assertTrue(fi.ctrl_bases["A"] == 2)
        self.assertTrue(fi.ctrl_bases["a"] == 0)
        
    def test13(self):
        fi = FisherInfo()
        fi.ref = "G"
        mp_list = []
        mp_list.append("chr1")
        mp_list.append("10000")
        mp_list.append("10001")
        mp_list.append("")
        mp_list.append("$.....AAAA,,a")
        mp_list.append("IIIIIIIIIIII")
        mp_list.append("")
        mp_list.append("......AA")
        mp_list.append("IIIIIIII")
        mp_list.append("")
        mp_list.append("AAAAAaaaaTT")
        mp_list.append("IIIIIIIIIII")
        fi.set_mpileup_data(mp_list)
        self.assertTrue(fi.tumor_bases["G"] == 5)
        self.assertTrue(fi.tumor_bases["g"] == 2)
        self.assertTrue(fi.tumor_bases["A"] == 4)
        self.assertTrue(fi.tumor_bases["a"] == 1)
        self.assertTrue(fi.ctrl_bases["G"] == 6)
        self.assertTrue(fi.ctrl_bases["g"] == 0)
        self.assertTrue(fi.ctrl_bases["A"] == 2)
        self.assertTrue(fi.ctrl_bases["a"] == 0)
        self.assertTrue(fi.rna_bases["G"] == 0)
        self.assertTrue(fi.rna_bases["g"] == 0)
        self.assertTrue(fi.rna_bases["A"] == 5)
        self.assertTrue(fi.rna_bases["a"] == 4)
        self.assertTrue(fi.rna_bases["T"] == 2)
        
    def test21(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 0
        fi.tumor_bases["a"] = 1
        fi.ref = "G"
        ret = fi.get_tumor_misrate("A")
        self.assertTrue(ret == 0.1)

    def test22(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 0
        fi.tumor_bases["a"] = 1
        fi.ctrl_bases["G"] = 0
        fi.ctrl_bases["g"] = 0
        fi.ctrl_bases["A"] = 0
        fi.ctrl_bases["a"] = 0
        fi.ref = "G"
        ret = fi.get_ctrl_misrate("A")
        self.assertTrue(ret == 0.0)
        
    def test23(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 0
        fi.tumor_bases["a"] = 1
        fi.ctrl_bases["G"] = 0
        fi.ctrl_bases["g"] = 0
        fi.ctrl_bases["A"] = 0
        fi.ctrl_bases["a"] = 0
        fi.rna_bases["G"] = 1
        fi.rna_bases["g"] = 1
        fi.rna_bases["T"] = 1
        fi.rna_bases["t"] = 1
        fi.rna_bases["a"] = 1
        fi.ref = "G"
        ret = fi.get_rna_misrate("A")
        self.assertTrue(ret == 0.2)
        
    def test31(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 1
        fi.tumor_bases["a"] = 0
        fi.ref = "G"
        ret = fi.get_tumor_strand_ratio("A")
        self.assertTrue(ret == 1)
        
    def test32(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 1
        fi.tumor_bases["a"] = 1
        fi.ref = "G"
        ret = fi.get_tumor_strand_ratio("A")
        self.assertTrue(ret == 0.5)

    def test33(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 1
        fi.tumor_bases["a"] = 1
        fi.ref = "G"
        ret = fi.get_ctrl_strand_ratio("A")
        self.assertTrue(ret == -1)
        
    def test34(self):
        fi = FisherInfo()
        fi.tumor_bases["G"] = 5
        fi.tumor_bases["g"] = 4
        fi.tumor_bases["A"] = 1
        fi.tumor_bases["a"] = 1
        fi.rna_bases["A"] = 0
        fi.rna_bases["a"] = 1
        fi.ref = "G"
        ret = fi.get_rna_strand_ratio("A")
        self.assertTrue(ret == 0)
        

