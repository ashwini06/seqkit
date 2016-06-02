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
                       '#SBATCH -p core -n 3 \n'
                       '#SBATCH -t 4:00:00\n'
                       '#SBATCH -e '+proj_dir+'/{treat}/scripts/{name}_postqc.stderr\n'
                       '#SBATCH -o '+proj_dir+'/{treat}/scripts/{name}_postqc.stdout\n'
                       '#SBATCH --mail-type=FAIL\n'
                       '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                       'module load bioinfo-tools\n'
                       'module load deepTools/2.2.3\n\n'
                       #'bamCompare --bamfile1 {treatment} --bamfile2 {control} --binSize 25 --ratio log2 --scaleFactorsMethod SES -o {postqc_dir}/{treat}_Vs_{ctrl}_log2ratio.bw\n'
                       'bamCoverage --bam {treatment} --binSize 25 --normalizeUsingRPKM -o {postqc_dir}/{treat}_coverage.bw -bl /home/ashwini/mm10_blacklisted-regions.bed\n'
                       'bamCoverage --bam {control} --binSize 25 --normalizeUsingRPKM -o {postqc_dir}/{ctrl}_coverage.bw -bl /home/ashwini/mm10_blacklisted-regions.bed\n'
                       'plotFingerprint -b {treatment} {control} -plot {postqc_dir}/{treat}_Vs_{ctrl}_fingerprint.png --labels {treat} {ctrl}\n'
                       'multiBamSummary bins --bamfiles '+proj_dir+'/*/alignment_bowtie/bam_files/*_rmdup.bam -out {postqc_dir}/results.npz\n'
                       'plotCorrelation --corData {postqc_dir}/results.npz --plotFile {postqc_dir}/scatterplot.pdf --corMethod pearson --whatToPlot scatterplot --skipZeros\n'
                       'computeMatrix scale-regions -S {postqc_dir}/{treat}_coverage.bw {postqc_dir}/{ctrl}_coverage.bw -R {ucsc_file} --skipZeros -o {postqc_dir}/matrix.mat.gz --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000\n'
                       'plotHeatmap -m {postqc_dir}/matrix.mat.gz -out {postqc_dir}/{treat}_Vs_{ctrl}_heatmap.png\n')


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
        job_file = os.path.join(proj_dir,treat,"{}/{}_{}.sh".format("scripts",name,"postqc"))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(sbatch_template.format(treat=treat, ctrl=ctrl, name=name, treatment=treat_fl, control=control_fl, ucsc_file=ucsc_file, postqc_dir=postqc_dir))
#        subprocess.check_call(['sbatch',job_file])


