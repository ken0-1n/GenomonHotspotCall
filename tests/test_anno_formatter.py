
import sys
import unittest, subprocess
import os, tempfile, shutil, filecmp
from hotspot_call import anno_formatter as anf
from hotspot_call.fisher_info import FisherInfo


class TestANNOFormatter(unittest.TestCase):

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
        ret = anf.make_record(fi, alt, is_rna, is_ctrl)
        ans = "chr1\t10000\t10001\tC\tA\t200\t100\t41\t1\t120,60,80,40\t30,0,11,1\t100,100,0,0\t1,40,0,0\t0.500\t0.600\t0.024\t0.000\t9.021\t6.680"
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
        ret = anf.make_record(fi, alt, is_rna, is_ctrl)
        ans = "chr1\t10000\t10001\tC\tA\t200\t100\t41\t1\t120,60,80,40\t30,0,11,1\t100,100,0,0\t1,40,0,0\t0.500\t0.600\t0.024\t0.000\t9.021\t6.680\t39\t13\t30,10,9,3\t0.333"
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
        ret = anf.make_record(fi, alt, is_rna, is_ctrl)
        ans = "chr1\t10000\t10001\tC\tA\t200\t100\t41\t1\t120,60,80,40\t100,100,0,0\t0.500\t0.600\t0.455\t0.500\t0.545\t6.680"
        self.assertTrue(ret == ans)
   
