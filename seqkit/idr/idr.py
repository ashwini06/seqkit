#!/usr/bin/env python
""" Methods and functionalites to perform peak-calling """
import os
import subprocess
from glob import glob
from seqkit import CONFIG as conf


def run_dr(project, input_file):
    """ Will run the idr analysis to check for biological replicate consistence """
    root_dir = conf.get('root_dir','')  
    proj_dir = os.path.join (root_dir,project)
    sbatch_template = ('#!/bin/bash -l\n'
                       '#SBATCH -A b2012025\n'
                       '#SBATCH -J {name}_peakcall\n'
                       '#SBATCH -p core -n 1 \n'
                       '#SBATCH -t 5:00:00\n'
                       '#SBATCH -o '+proj_dir+'/{rep1}_Vs_{rep2}/scripts/{name}_idr.stdout\n'
                       '#SBATCH -e '+proj_dir+'/{rep1}_Vs_{rep2}/scripts/{name}_idr.stderr\n'
                       '#SBATCH --mail-type=FAIL\n'
                       '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                       'module load bioinfo-tools\n'
                       'sort -k 8,8nr {rep1_dir}/*.narrowPeak > {rep1_dir}/tmp.regionPeak\n'
                       'intersectBed -a {rep1_dir}/tmp.regionPeak -b {mm10_blacklisted-regions.bed} > {rep1_dir}/cleanedpeaks.regionPeak\n'
                       'sort -k 8,8nr {rep2_dir}/*.narrowPeak > {rep2_dir}/tmp.regionPeak\n'
                       'intersectBed -a {rep2_dir}/tmp.regionPeak -b {mm10_blacklisted-regions.bed} > {rep2_dir}/cleanedpeaks.regionPeak\n' 
                        'Rscript batch-consistency-analysis.r rep1,rep2 -1 idr_op 0 F p.value {mm10.genome}\n'

)


    pk_file = open(input_file,'r')
    pk_file.next()
    for ln in iter(pk_file):    
        ln = ln.strip()
        ln =  ln.split('\t')
        rep1 = ln[0]
        rep2 = ln[1]
        name = "{}_Vs_{}".format(rep1,rep2)
        rep1_dir = ''.join(glob("{}/{}/macs2_*".format(proj_dir,rep1)))
        peaks_dir = os.path.join(proj_dir,treat,"{}_{}".format(peak_call,mode))
        if not os.path.exists(peaks_dir):
            os.makedirs(peaks_dir)
        job_fl = os.path.join(proj_dir,treat,"scripts","{}_peakcall.sh".format(name))
        template_pc = sbatch_template + template
        with open(job_fl,'w') as jb_fl:
            jb_fl.write(template_pc.format(name=name,treat=treat, treatment=treat_fl, control=control_fl,peaks_dir=peaks_dir))
    subprocess.check_call(['sbatch',job_fl]) 
