#!/usr/bin/env python
""" Methods and functionalites to perform peak-calling """
import os
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit import utils

def run_denovo(project,peak_call,slurm=False,job_file=None):
    """ Will run the de-novo motif analysis """
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir,project)
    sample = map(str,glob(os.path.join(proj_dir, "*", "*{}*".format(peak_call),"*xls")))
    motif_r = os.path.join(os.path.dirname(utils.__file__),"motifanalysis.r")
    sbatch_template = ('#!/bin/bash -l\n'
    '#SBATCH -A b2012025\n'
    '#SBATCH -J {nm}_motifanalysis\n'
    '#SBATCH -p core -n 1 \n'
    '#SBATCH -t 3:00:00\n'
    '#SBATCH --mail-type=FAIL\n'
    '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n')
    template_denovo = ('\n## Running de-nove motif analysis\n'
    'module load bioinfo-tools\n'
    'module load MEMEsuite/4.11.1\n'    
    'Rscript '+motif_r+' {ip_fl} {op_dir} {op_fl}\n'
    )

    for xls in sample:
        nm = os.path.basename(xls).replace(".xls","")
        op_dir = os.path.join(os.path.dirname(os.path.dirname(xls)),"motif")
        if not os.path.exists(op_dir):
            os.mkdir(op_dir)
        op_fl = nm +("_seq.fa")
        job_file = os.path.join(os.path.dirname(os.path.dirname(xls)),"scripts","{}_denovo.sh".format(nm))
        template = sbatch_template+template_denovo
        with open(job_file,'w') as jb_fl:
            jb_fl.write(template.format(ip_fl=xls,op_dir=op_dir,op_fl=op_fl,nm=nm))




