#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
import pdb
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def bamcov(project, genefile, sample = None):
    """Will run the postqc"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    sbatch_template = ('#!/bin/bash -l\n'
                       '#SBATCH -A b2012025\n'
                       '#SBATCH -J {sam}_postqc\n'
                       '#SBATCH -p core -n 1 \n'
                       '#SBATCH -t 4:00:00\n'
                       '#SBATCH -e {sam_dir}/scripts/{sam}_postqc.stderr\n'
                       '#SBATCH -o {sam_dir}/scripts/{sam}_postqc.stdout\n'
                       '#SBATCH --mail-type=FAIL\n'
                       '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                       'module load bioinfo-tools\n'
                       'module load deepTools/2.0.1\n'
                       'bamCoverage -b {ipfile} -o {outfile} -normalizeUsingRPKM=True --extendReads\n'
                       #'bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw\n'
                       'computeMatrix scale-regions -S {outfile} -R {ucsc_file} --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --skipZeros -o {matrix_fl} --outFileSortedRegions {sorted_fl}\n'
                       'plotHeatmap -m {matrix_fl} -out {hmap}\n')

    if sample:
        if os.path.isdir(os.path.join(proj_dir, sample)):
            samples = [sample]
        else:
            raise SystemExit("Given sample {} is not found in project directory {}".format(sample, proj_dir))
    else:
        samples = find_samples(proj_dir)
    
    ucsc_file = genefile

    for sam in samples:
        sam_dir = os.path.join(proj_dir, sam)
        src_dir = os.path.join(sam_dir, 'scripts')
        if not os.path.exists(src_dir):
            os.mkdir(src_dir)
        postqc_dir = os.path.join(sam_dir,"deepTools")
        if not os.path.exists(postqc_dir):
            os.makedirs(postqc_dir)
        ip_bam = ''.join(glob("{}/{}/{}/*bam".format(sam_dir,"alignment_*","bam_files")))
        name=os.path.basename(ip_bam).replace("_sorted.bam","")
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
        job_file = os.path.join(src_dir, "{}_{}.sh".format(name,"postqc"))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(sbatch_template.format(sam=sam, sam_dir=sam_dir, postqc_dir=postqc_dir, ipfile=ip_bam, outfile=op_fl,matrix_fl=matrix_fl, sorted_fl=sorted_fl,hmap=hmap, ucsc_file=ucsc_file))
        #subprocess.check_call(['sbatch',job_file])


