#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
import pdb
import subprocess
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def run_qc(project):
    """Will run the QC"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    fastqc_sbatch_template = ('#!/bin/bash -l\n'
                              '#SBATCH -A b2012025\n'
                              '#SBATCH -J {sam}_fastqc\n'
                              '#SBATCH -p core -n 2 \n' 
                              '#SBATCH -t 2:00:00\n'
                              '#SBATCH -e {sam_dir}/scripts/{sam}_fastqc.stderr\n'
                              '#SBATCH -o {sam_dir}/scripts/{sam}_fastqc.stdout\n'
                              '#SBATCH --mail-type=FAIL\n'
                              '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                              'module load bioinfo-tools\n'
                              'module load FastQC/0.11.5\n'
                              'cd '+proj_dir+'\n'
                              'fastqc -o {fastqc_dir} -f fastq {fq_files}\n')
    samples = find_samples(proj_dir,file_type="fastq")
    for sam in samples.keys():
        fq_fls = samples[sam]
        pdb.set_trace()
        sam_dir = os.path.join(proj_dir, sam)
        src_dir = os.path.join(sam_dir, 'scripts')
        if not os.path.exists(src_dir):
            os.mkdir(src_dir)
        fastqc_dir = os.path.join(sam_dir,'fastqc')
        if not os.path.exists(fastqc_dir):
            os.mkdir(fastqc_dir)
        job_file = os.path.join(src_dir, "{}_fastqc.sh".format(sam))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(fastqc_sbatch_template.format(sam=sam, sam_dir=sam_dir,fastqc_dir=fastqc_dir, fq_files=" ".join(fq_fls)))
        subprocess.check_call(['sbatch',job_file])

