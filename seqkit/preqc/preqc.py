#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os

def run_qc(output=None):
    """Will run the QC"""
    if not output:
        output = os.getcwd()
    print "Under development, when done with run the QC and save in {}".format(output)
