#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
import pdb
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def bamcov(project, genefile, input_file):
    """Will run the postqc"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    sbatch_template = ('#!/bin/bash -l\n'
                       '#SBATCH -A b2012025\n'
                       '#SBATCH -J {name}_postqc\n'
                       '#SBATCH -p core -n 1 \n'
                       '#SBATCH -t 4:00:00\n'
                       '#SBATCH -e '+proj_dir+'/{treat}/scripts/{name}_postqc.stderr\n'
                       '#SBATCH -o '+proj_dir+'/{treat}/scripts/{name}_postqc.stdout\n'
                       '#SBATCH --mail-type=FAIL\n'
                       '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                       'module load bioinfo-tools\n'
                       'module load deepTools/2.0.1\n'
                       'bamCompare -b1 {treatment} -b2 {control} -o {outfile}\n'
                       'computeMatrix scale-regions -S {outfile} -R {ucsc_file} --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o {matrix_fl} --outFileSortedRegions {sorted_fl}\n'
                       #plotting heatmap
                        'plotHeatmap -m {matrix_fl} -out {hmap}\n')

    ucsc_file = genefile

    pk_file = open(input_file,'r')
    pk_file.next()
    for ln in iter(pk_file):
        ln = ln.strip()
        ln =  ln.split('\t')
        treat = ln[0]
        ctrl = ln[1]
        name = "{}_Vs_{}".format(treat,ctrl)
        treat_fl = ''.join(glob("{}/{}/alignment_*/bam_files/{}*sorted_rmdup.bam".format(proj_dir,treat,treat)))
        control_fl = ''.join(glob("{}/{}/alignment_*/bam_files/{}*sorted_rmdup.bam".format(proj_dir,ctrl,ctrl)))
        postqc_dir = os.path.join(proj_dir,treat,"deepTools")
        if not os.path.exists(postqc_dir):
            os.mkdir(postqc_dir)
        op_fl = os.path.join(postqc_dir,"{}_{}.bw".format(name,"coverage"))
        if not os.path.exists(op_fl):
            os.mknod(op_fl)    
        matrix_fl = os.path.join(postqc_dir,"{}_{}.gz".format(name,"mat"))
        if not os.path.exists(matrix_fl):
            os.mknod(matrix_fl)
        sorted_fl = os.path.join(postqc_dir,"{}_{}.bed".format(name,"sorted_genes"))
        if not os.path.exists(sorted_fl):
            os.mknod(sorted_fl)
        hmap = os.path.join(postqc_dir,"{}_{}.png".format(name,"heatmap"))
        if not os.path.exists(hmap):
            os.mknod(hmap)
        job_file = os.path.join(proj_dir,treat,"{}/{}_{}.sh".format("scripts",name,"postqc"))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(sbatch_template.format(treat=treat,name=name,treatment=treat_fl, control=control_fl,outfile=op_fl,matrix_fl=matrix_fl, sorted_fl=sorted_fl,hmap=hmap, ucsc_file=ucsc_file))
        subprocess.check_call(['sbatch',job_file])


