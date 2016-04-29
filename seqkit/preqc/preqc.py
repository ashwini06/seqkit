#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
import pdb
from seqkit.utils.find_samples import find_samples

def run_qc(project):
    """Will run the QC"""
    samples = find_samples(project)
	for sam in samples.keys():
		fq_fls = samples[sam]
		pdb.set_trace()
    	## do everything after here ##
		fl = open(+sam,'_fastqc.sh','w')
		fl.write(
		'#!/bin/bash -l\n'
		'#SBATCH -A b2012025\n'
		'#SBATCH -J '+sam''_fastqc'\n'
		'#SBATCH -p node -N 1 \n'	
		'#SBATCH -t 2:00:00\n'
		'#SBATCH --mail-type=FAIL\n'
		'#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
		'module load bioinfo-tools\n'
		'module load FastQC/0.11.5\n'
		#'op_dir=\'\'\n'
		#'fastqc -o '+outdir' -f fastq '+fq_fls'\n'
		'fastqc -f fastq '+fq_fls'\n'
		)
	fl.close()	
	subprocess.check_call(['sbatch', +sam,'_fastqc.sh'])