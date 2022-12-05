#! /usr/bin/env python
import sys, os
from .fisher_info import FisherInfo


def make_record(fi, alt, is_rna, is_ctrl):

    outstr = (
    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
    fi.chr,
    fi.start,
    fi.end,
    fi.ref,
    alt,
    fi.get_tumor_depth(),
    fi.get_tumor_base_total(alt),
    fi.get_ctrl_depth(),
    fi.get_ctrl_base_total(alt))
    + '\t{0},{1},{2},{3}'.format(
    fi.get_tumor_depth_plus_strand(),
    fi.get_tumor_base_plus_strand(alt),
    fi.get_tumor_depth_minus_strand(),
    fi.get_tumor_base_minus_strand(alt))
    )

    if is_ctrl:
        outstr = outstr + '\t{0},{1},{2},{3}'.format(
        fi.get_ctrl_depth_plus_strand(),
        fi.get_ctrl_base_plus_strand(alt),
        fi.get_ctrl_depth_minus_strand(),
        fi.get_ctrl_base_minus_strand(alt))

    outstr = outstr + '\t{0},{1},{2},{3}'.format(
    fi.tumor_bases["A"] + fi.tumor_bases["a"],
    fi.tumor_bases["C"] + fi.tumor_bases["c"],
    fi.tumor_bases["G"] + fi.tumor_bases["g"],
    fi.tumor_bases["T"] + fi.tumor_bases["t"])

    if is_ctrl:
        outstr = outstr + '\t{0},{1},{2},{3}'.format(
        fi.ctrl_bases["A"] + fi.ctrl_bases["a"],
        fi.ctrl_bases["C"] + fi.ctrl_bases["c"],
        fi.ctrl_bases["G"] + fi.ctrl_bases["g"],
        fi.ctrl_bases["T"] + fi.ctrl_bases["t"])

    outstr = outstr + '\t{0:.3f}'.format(fi.get_tumor_misrate(alt))
    outstr = outstr + '\t{0:.3f}'.format(fi.get_tumor_strand_ratio(alt))

    if is_ctrl:
        outstr = outstr + '\t{0:.3f}'.format(fi.get_ctrl_misrate(alt))
        ctrl_strand_ratio = fi.get_ctrl_strand_ratio(alt)
        if ctrl_strand_ratio > 0:
            outstr = outstr + '\t{0:.3f}'.format(ctrl_strand_ratio)
        else:
            outstr = outstr + "\t---"
        outstr = outstr + '\t{0:.3f}'.format(fi.get_fisher_pvalue(alt))

    outstr = outstr + '\t{0:.3f}'.format(fi.get_lod_score(alt))

    if is_rna:
        outstr = outstr + (
        '\t{0}\t{1}'.format(
        fi.get_rna_depth(),
        fi.get_rna_base_total(alt))
        + '\t{0},{1},{2},{3}'.format(
        fi.get_rna_depth_plus_strand(),
        fi.get_rna_base_plus_strand(alt),
        fi.get_rna_depth_minus_strand(),
        fi.get_rna_base_minus_strand(alt))
        )
        outstr = outstr + '\t{0:.3f}'.format(fi.get_rna_misrate(alt))
    return outstr

