#!/usr/bin/env python

import os
from glob import glob

def find_samples(proj_dir, file_type="fastq", aligner="bowtie2"):
    """Find the samples in a given proj dir"""
    current_dir = os.getcwd()
    samples_fls = {}
    os.chdir(proj_dir)
    samples = [sam for sam in os.listdir(os.getcwd()) if os.path.isdir(sam)]
    for sam in samples:
        if file_type="fastq"
            samples_fls[sam] = glob("{}/Rawdata/*.fastq*".format(sam))
        elif file_type="bam"
            samples_fls[sam] = glob("{}/alignment_{}/bam_files/*.bam".format(sam, aligner))
    os.chdir(current_dir)
    return samples_fls
