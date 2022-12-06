
import sys
import unittest, subprocess
import os, tempfile, shutil, filecmp
from hotspot_call import process_mutation as pm
from hotspot_call.fisher_info import FisherInfo


class TestProcessMutation(unittest.TestCase):

    # get_ebpval():
    def test01(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp() 

        hotspot_file = tmp_dir+"/hotspot.txt"
        with open(hotspot_file, "w") as hout:
            print('chr1\t100\t100\tA\tC\tMTOR\tc.10T>G\tp.I30M\tSNP', file=hout)

        d_hotspot = pm.read_hotspot_file(hotspot_file)
        self.assertTrue(d_hotspot['chr1\t100\t100\tA'] == "C")

    def test02(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp() 

        hotspot_file = tmp_dir+"/hotspot.txt"
        with open(hotspot_file, "w") as hout:
            print('chr1\t100\t100\tAA\tC\tMTOR\tc.10T>G\tp.I30M\tSNP', file=hout)

        d_hotspot = pm.read_hotspot_file(hotspot_file)
        self.assertTrue(len(d_hotspot) == 0)

    def test03(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp() 

        hotspot_file = tmp_dir+"/hotspot.txt"
        with open(hotspot_file, "w") as hout:
            print('chr1\t100\t100\tA\tCR\tMTOR\tc.10T>G\tp.I30M\tSNP', file=hout)

        d_hotspot = pm.read_hotspot_file(hotspot_file)
        self.assertTrue(len(d_hotspot) == 0)

    def test11(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        test_file = cur_dir + "/../data/test11.txt"
        answer_file = cur_dir + "/../data/ans11.txt"
        with open(test_file, "w") as hout:
            pm.print_anno_header(False, hout, True) 
        self.assertTrue(filecmp.cmp(test_file, answer_file, shallow=False))

    def test12(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        test_file = cur_dir + "/../data/test12.txt"
        answer_file = cur_dir + "/../data/ans12.txt"
        with open(test_file, "w") as hout:
            pm.print_anno_header(False, hout, False) 
        self.assertTrue(filecmp.cmp(test_file, answer_file, shallow=False))

    def test13(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        test_file = cur_dir + "/../data/test13.txt"
        answer_file = cur_dir + "/../data/ans13.txt"
        with open(test_file, "w") as hout:
            pm.print_anno_header(True, hout, True) 
        self.assertTrue(filecmp.cmp(test_file, answer_file, shallow=False))

    def test21(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_fa = os.environ['HOME'] + "/database/v0/Homo_sapiens_assembly38.fasta"
        ref_dict = os.environ['HOME'] + "/database/v0/Homo_sapiens_assembly38.dict"
        test_file = cur_dir + "/../data/test21.txt"
        answer_file = cur_dir + "/../data/ans21.txt"
        with open(test_file, "w") as hout:
            pm.print_vcf_header(False, "Sample1", "Sample2", None, ref_fa, ref_dict, hout, True)
        self.assertTrue(filecmp.cmp(test_file, answer_file, shallow=False))

    def test22(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_fa = os.environ['HOME'] + "/database/v0/Homo_sapiens_assembly38.fasta"
        ref_dict = os.environ['HOME'] + "/database/v0/Homo_sapiens_assembly38.dict"
        test_file = cur_dir + "/../data/test22.txt"
        answer_file = cur_dir + "/../data/ans22.txt"
        with open(test_file, "w") as hout:
            pm.print_vcf_header(False, "Sample1", "Sample2", None, ref_fa, ref_dict, hout, False)
        self.assertTrue(filecmp.cmp(test_file, answer_file, shallow=False))

    def test23(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))

        ref_fa = os.environ['HOME'] + "/database/v0/Homo_sapiens_assembly38.fasta"
        ref_dict = os.environ['HOME'] + "/database/v0/Homo_sapiens_assembly38.dict"
        test_file = cur_dir + "/../data/test23.txt"
        answer_file = cur_dir + "/../data/ans23.txt"
        with open(test_file, "w") as hout:
            pm.print_vcf_header(True, "Sample1", "Sample2", "Sample3", ref_fa, ref_dict, hout, True)
        self.assertTrue(filecmp.cmp(test_file, answer_file, shallow=False))

    def test31(self):
        mutReg = "chr1:100-200"
        m_params = "-test option".split(" ")
        is_rna = False
        is_ctrl = True
        bam_tumor = "data/tumor.cram"
        bam_control = "data/ctrl.cram"
        bam_rna = None
        command = pm.make_pileup_cmmand(mutReg, m_params, is_rna, is_ctrl, bam_tumor, bam_control, bam_rna)
        ans = ["samtools", "mpileup", "-r", "chr1:100-200", "-test", "option", "data/tumor.cram", "data/ctrl.cram"]
        self.assertTrue(command == ans)

    def test32(self):
        mutReg = "chr1:100-200"
        m_params = "-test option".split(" ")
        is_rna = False
        is_ctrl = False
        bam_tumor = "data/tumor.cram"
        bam_control = None
        bam_rna = None
        command = pm.make_pileup_cmmand(mutReg, m_params, is_rna, is_ctrl, bam_tumor, bam_control, bam_rna)
        ans = ["samtools", "mpileup", "-r", "chr1:100-200", "-test", "option", "data/tumor.cram"]
        self.assertTrue(command == ans)

    def test33(self):
        mutReg = "chr1:100-200"
        m_params = "-test option".split(" ")
        is_rna = True
        is_ctrl = True
        bam_tumor = "data/tumor.cram"
        bam_control = "data/ctrl.cram"
        bam_rna = "data/rna.cram"
        command = pm.make_pileup_cmmand(mutReg, m_params, is_rna, is_ctrl, bam_tumor, bam_control, bam_rna)
        ans = ["samtools", "mpileup", "-r", "chr1:100-200", "-test", "option", "data/tumor.cram", "data/ctrl.cram", "data/rna.cram"]
        self.assertTrue(command == ans)

    def test41(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 6
        fi.ctrl_bases["c"] = 4
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.5
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = True
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == False)
   
    def test42(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 6
        fi.ctrl_bases["c"] = 4
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.51
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = True
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == True)
   
    def test43(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 6
        fi.ctrl_bases["c"] = 4
        fi.ctrl_bases["a"] = 1
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.50
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = True
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == True)
   
    def test44(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 6
        fi.ctrl_bases["c"] = 4
        alt = "A"
        min_lod_score = 6.7
        min_tumor_misrate = 0.50
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = True
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == True)

    def test45(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 6
        fi.ctrl_bases["c"] = 4
        fi.ctrl_bases["a"] = 1
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.50
        max_ctrl_misrate = 0.1
        ratio_ctrl = 0.1
        is_ctrl = True
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == True)
      
    def test46(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 0
        fi.ctrl_bases["c"] = 0
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.5
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = False
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == False)
   
    def test47(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 0
        fi.ctrl_bases["c"] = 0
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.51
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = False
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == True)
   
    def test48(self):
        fi = FisherInfo()
        fi.tumor_quals["A"] = [10,10,10,10,10,10,10]
        fi.tumor_bases["A"] = 6
        fi.tumor_bases["a"] = 4 
        fi.tumor_bases["C"] = 6
        fi.tumor_bases["c"] = 4
        fi.ctrl_bases["C"] = 0
        fi.ctrl_bases["c"] = 0
        fi.ctrl_bases["a"] = 1
        alt = "A"
        min_lod_score = 6.6
        min_tumor_misrate = 0.50
        max_ctrl_misrate = 0.0
        ratio_ctrl = 0.1
        is_ctrl = False
        ret = pm.is_filter_fisher_info(fi, alt, min_lod_score, min_tumor_misrate, max_ctrl_misrate, ratio_ctrl, is_ctrl)
        self.assertTrue(ret == False)
   

    # if fi.get_lod_score(alt) < min_lod_score or fi.get_tumor_misrate(alt) < min_tumor_misrate:
    # if fi.get_ctrl_misrate(alt) > max_ctrl_misrate or fi.get_ctrl_misrate(alt) > (fi.get_tumor_misrate(alt) * ratio_ctrl):

