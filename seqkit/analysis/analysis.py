#!/usr/bin/env python
"""Methods and functionalites to do analysis"""
import os
import pdb
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def run_bowtie(project):
    """Will run the bowtie-alignment"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    bowtie2_index='/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index'
    bowtie2_template = ('#!/bin/bash -l\n'
    		           '#SBATCH -A b2012025\n'
    		           '#SBATCH -J {sam}_fastqc\n'
    		           '#SBATCH -p core -n 2 \n'	
    		           '#SBATCH -t 2:00:00\n'
    		           '#SBATCH --mail-type=FAIL\n'
    		           '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
    		           'module load bioinfo-tools\n'
    		           'module load bowtie2/2.2.6\n'
    		           'module load samtools/1.3\n'
		 	   'cd '+proj_dir+'\n'
			   'for run in {fq_files};do\n'
			   'nm=$(basename \$run)\n'
			   'nm=\${nm/.fastq/}\n'
			   'bowtie2 -t -p 8 '+bowtie_index+' -q ${run} -S {bowtie_dir}/${nm}.sam > {bowtie_dir}/${nm}_bowtie.log\n'
			   'done\n'
                        )
    samples = find_samples(proj_dir)
    for sam in samples.keys():
        fq_fls = samples[sam]
	#pdb.set_trace()
        sam_dir = os.path.join(proj_dir, sam)
	src_dir = os.path.join(sam_dir, 'scripts')
        if not os.path.exists(src_dir):
        	os.mkdir(src_dir)
	bowtie2_dir = os.path.join(sam_dir,'bowtie2')
	if not os.path.exists(bowtie2_dir):
		os.mkdir(bowtie2_dir)
        job_file = os.path.join(src_dir, "{}_bowtie2.sh".format(sam))
	with open(job_file, 'w') as jb_fl:
        	jb_fl.write(bowtie2_template.format(sam=sam, sam_dir=sam_dir,bowtie2_dir=bowtie2_dir, fq_files=" ".join(fq_fls)))
	#subprocess.check_call(['sbatch',job_file])
