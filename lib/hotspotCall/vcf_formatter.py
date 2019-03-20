#! /usr/bin/env python
import sys, os
from fisher_info import FisherInfo

def make_record(fi, alt, is_rna):
    
    outstr = (
    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\tFP={7:.3f};LS={8:.3f}\tDP:DPF:DPR:AD:ADF:ADR:AF:SB'.format(
    fi.chr,
    fi.start,
    '.',
    fi.ref,
    alt,
    '.',
    '.',
    fi.get_fisher_pvalue(alt),
    fi.get_lod_score(alt)
    )
    + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6:.3f}:{7:.3f}'.format(
    fi.get_tumor_depth(),
    fi.get_tumor_depth_plus_strand(),
    fi.get_tumor_depth_minus_strand(),
    fi.get_tumor_base_total(alt),
    fi.get_tumor_base_plus_strand(alt),
    fi.get_tumor_base_minus_strand(alt),
    fi.get_tumor_misrate(alt),
    fi.get_tumor_strand_ratio(alt),
    )
    + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6:.3f}'.format(
    fi.get_ctrl_depth(),
    fi.get_ctrl_depth_plus_strand(),
    fi.get_ctrl_depth_minus_strand(),
    fi.get_ctrl_base_total(alt),
    fi.get_ctrl_base_plus_strand(alt),
    fi.get_ctrl_base_minus_strand(alt),
    fi.get_ctrl_misrate(alt)
    ))
    ctrl_strand_ratio = fi.get_ctrl_strand_ratio(alt)
    if ctrl_strand_ratio > 0:
        outstr = outstr + ':{0:.3f}'.format(ctrl_strand_ratio)
    else:
        outstr = outstr + ":."

    if is_rna:
        outstr = outstr + (
        '\t{0}:{1}:{2}:{3}:{4}:{5}:{6:.3f}'.format(
        fi.get_rna_depth(),
        fi.get_rna_depth_plus_strand(),
        fi.get_rna_depth_minus_strand(),
        fi.get_rna_base_total(alt),
        fi.get_rna_base_plus_strand(alt),
        fi.get_rna_base_minus_strand(alt),
        fi.get_rna_misrate(alt)
        ))
        rna_strand_ratio = fi.get_rna_strand_ratio(alt)
        if rna_strand_ratio > 0:
            outstr = outstr + ':{0:.3f}'.format(rna_strand_ratio)
        else:
            outstr = outstr + ":."

    return outstr

