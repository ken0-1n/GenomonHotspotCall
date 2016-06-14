#! /usr/bin/env python
import sys, os, re
from fisher_info import FisherInfo

def make_record(fi, alt):

    outstr = (
    '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}'.format(
    fi.chr,
    fi.start,
    fi.end,
    fi.ref,
    alt,
    fi.get_target_depth_total(),
    fi.get_target_total(alt),
    fi.get_control_depth_total(),
    fi.get_control_total(alt))
    + '\t{0},{1},{2},{3}'.format(
    fi.get_target_depth_plus_strand(),
    fi.get_target_plus_strand(alt),
    fi.get_target_depth_minus_strand(),
    fi.get_target_minus_strand(alt))
    + '\t{0},{1},{2},{3}'.format(
    fi.get_control_depth_plus_strand(),
    fi.get_control_plus_strand(alt),
    fi.get_control_depth_minus_strand(),
    fi.get_control_minus_strand(alt))
    + '\t{0},{1},{2},{3}'.format(
    fi.get_target_A_total(),
    fi.get_target_C_total(),
    fi.get_target_G_total(),
    fi.get_target_T_total())
    + '\t{0},{1},{2},{3}'.format(
    fi.get_control_A_total(),
    fi.get_control_C_total(),
    fi.get_control_G_total(),
    fi.get_control_T_total())
    )
    outstr = outstr + '\t{0:.3f}'.format(fi.get_target_misrate(alt))
    t_ratio = fi.get_target_strand_ratio(alt)
    if t_ratio >= 0:
       outstr = outstr + '\t{0:.3f}'.format(t_ratio)
    else:
       outstr = outstr + "\t---"
    outstr = outstr + '\t{0:.3f}'.format(fi.get_control_misrate(alt))
    t_ratio = fi.get_control_strand_ratio(alt)
    if t_ratio >= 0:
       outstr = outstr + '\t{0:.3f}'.format(t_ratio)
    else:
       outstr = outstr + "\t---"
    outstr = outstr + '\t{0:.3f}'.format(fi.get_fisher_pvalue(alt)) 
    
    return outstr

