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
                       #'bamCompare -b1 {treatment} -b2 {control} -o {outfile}\n\n'
                       'bamCompare --bamfile1 {treatment} --bamfile2 {control} --binSize 25 --ratio log2 --scaleFactorsMethod SES -o {postqc_dir}/{treat}_Vs_{ctrl}_log2ratio.bw\n'
                      #'computeMatrix scale-regions --regionsFileName {ucsc_file} --scoreFileName {outfile} --upstream 1000 --downstream 1000 --regionBodyLength 1000 --binSize 1 --sortRegions no --sortUsing mean --averageTypeBins mean --outFileName {matrix_fl}\n\n'
                       ##'computeMatrix scale-regions -S {outfile} -R {ucsc_file} --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o {matrix_fl} --outFileSortedRegions {sorted_fl}\n'
                       #plotting heatmap
                       # 'plotHeatmap -m {matrix_fl} -out {hmap}\n')
                       'plotFingerprint -b {treatment} {control} -plot {postqc_dir}/{treat}_Vs_{ctrl}_fingerprint.png --labels {treat} {ctrl}\n')
                       # multiBamSummary bins --bamfiles /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/*/alignment_bowtie/bam_files/*_rmdup.bam -out /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/results.npz
#                        plotCorrelation --corData /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/results.npz --plotFile /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/scatterplot.pdf --corMethod pearson --whatToPlot scatterplot --skipZeros
#                        computeMatrix scale-regions -S /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/treat.bw -R /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/chr1_ucsc.bed --skipZeros -o /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/matrix.mat.gz
#plotHeatmap -m /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/matrix.mat.gz -out /proj/b2012025/RAW_DATA/ChIP_Data/Ascl1_US/Mark_Mash1_s4/deepTools/heatmap1_chr.png


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
        #hmap = os.path.join(postqc_dir,"{}_{}.png".format(name,"heatmap"))
        #if not os.path.exists(hmap):
        #    os.mknod(hmap)
        job_file = os.path.join(proj_dir,treat,"{}/{}_{}.sh".format("scripts",name,"postqc"))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(sbatch_template.format(treat=treat,ctrl=ctrl,name=name,treatment=treat_fl, control=control_fl,outfile=op_fl,ucsc_file=ucsc_file,postqc_dir=postqc_dir))
#        subprocess.check_call(['sbatch',job_file])


