#! /usr/bin/env python
import sys, os
from .fisher_info import FisherInfo

def make_record(fi, alt, is_rna, is_ctrl):
    
    outstr = (
    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(
    fi.chr,
    fi.start,
    '.',
    fi.ref,
    alt,
    '.',
    '.'
    ))

    if is_ctrl:
        outstr = outstr + '\tFP={0:.3f}'.format(fi.get_fisher_pvalue(alt))
    else:
        b10, bm, b90 = fi.calc_btdtri(alt)
        outstr = outstr + '\tB10={0:.3f};BM={1:.3f};B90={2:.3f}'.format(b10, bm, b90)
        
    outstr = outstr + ';LS={0:.3f}'.format(fi.get_lod_score(alt))
    outstr = outstr + '\tDP:DPF:DPR:AD:ADF:ADR:AF:SB'

    outstr = outstr + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6:.3f}:{7:.3f}'.format(
    fi.get_tumor_depth(),
    fi.get_tumor_depth_plus_strand(),
    fi.get_tumor_depth_minus_strand(),
    fi.get_tumor_base_total(alt),
    fi.get_tumor_base_plus_strand(alt),
    fi.get_tumor_base_minus_strand(alt),
    fi.get_tumor_misrate(alt),
    fi.get_tumor_strand_ratio(alt),
    )

    if is_ctrl:
        outstr = outstr + '\t{0}:{1}:{2}:{3}:{4}:{5}:{6:.3f}'.format(
        fi.get_ctrl_depth(),
        fi.get_ctrl_depth_plus_strand(),
        fi.get_ctrl_depth_minus_strand(),
        fi.get_ctrl_base_total(alt),
        fi.get_ctrl_base_plus_strand(alt),
        fi.get_ctrl_base_minus_strand(alt),
        fi.get_ctrl_misrate(alt)
        )
        if fi.get_ctrl_misrate(alt) > 0:
            outstr = outstr + ':{0:.3f}'.format(fi.get_ctrl_strand_ratio(alt))
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
        if fi.get_rna_misrate(alt) > 0:
            outstr = outstr + ':{0:.3f}'.format(fi.get_rna_strand_ratio(alt))
        else:
            outstr = outstr + ":."

    return outstr

