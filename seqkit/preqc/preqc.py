#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
import pdb
from seqkit.utils.find_samples import find_samples

def run_qc(project):
    """Will run the QC"""
    samples = find_samples(project)
    pdb.set_trace()
    ## do everything after here ##
    cmd = ['module','load','FastQC/0.11.5']
	subprocess.check_all(cmd)
	cmd = ['fastqc', samples, '-o' ]
	subprocess.check_all(cmd)
