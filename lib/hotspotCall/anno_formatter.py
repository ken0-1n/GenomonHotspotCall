#! /usr/bin/env python
import sys, os
from fisher_info import FisherInfo


def make_record(fi, alt, is_rna):

    outstr = (
    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(
    fi.chr,
    fi.start,
    fi.end,
    fi.ref,
    alt,
    fi.get_tumor_depth(),
    fi.get_tumor_base_total(alt))
    + '\t{0},{1},{2},{3}'.format(
    fi.get_tumor_depth_plus_strand(),
    fi.get_tumor_base_plus_strand(alt),
    fi.get_tumor_depth_minus_strand(),
    fi.get_tumor_base_minus_strand(alt))
    + '\t{0},{1},{2},{3}'.format(
    fi.get_ctrl_depth_plus_strand(),
    fi.get_ctrl_base_plus_strand(alt),
    fi.get_ctrl_depth_minus_strand(),
    fi.get_ctrl_base_minus_strand(alt))
    )
    outstr = outstr + '\t{0:.3f}'.format(fi.get_tumor_misrate(alt))
    outstr = outstr + '\t{0:.3f}'.format(fi.get_ctrl_misrate(alt))
    outstr = outstr + '\t{0:.3f}'.format(fi.get_lod_score(alt))
    outstr = outstr + '\t{0:.3f}'.format(fi.get_score_median(alt))

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

