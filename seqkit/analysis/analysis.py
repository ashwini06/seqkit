#!/usr/bin/env python
"""Methods and functionalites to do analysis"""
import os

def run_analysis(output=None):
    """Will run the analysis"""
    if not output:
        output = os.getcwd()
    print "Under development, when done with run the pipeline and save in {}".format(output)
