#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
from seqkit.utils.find_samples import find_samples

def run_qc(project):
    """Will run the QC"""
    samples = find_samples(project)
    ## do everything after here ##
