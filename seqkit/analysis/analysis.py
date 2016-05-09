#!/usr/bin/env python
"""Methods and functionalites to do analysis"""
import os
import pdb
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def run_align(project,aligner):
	"""Will run the preferred-alignment"""
	root_dir = conf.get('root_dir','')
	proj_dir = os.path.join(root_dir, project)
	if aligner == "bwa":
		align_module = 'module load bwa/0.7.12\n'
		align_index = '/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa'
    		align_block = ('bwa aln {align_index} ${{run}} > {align_dir}/${{nm}}.sai\n'
                       		'bwa samse {align_index} {align_dir}/${{nm}}.sai ${{run}}.fastq | samtools view -Sb - > {align_dir}/${{nm}}.bam\n'
				'rm ${{nm}}.sai\n')
	elif aligner == "bowtie2":
		align_module = 'module load bowtie2/2.2.6\n'
		align_index = "/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome"
		align_block =  ('bowtie2 -t -p 8 -k2 --very-sensitive -x {align_index} -q ${{run}} -S {align_dir}/${{nm}}.sam > {align_dir}/${{nm}}_bowtie2.log\n'
                       'samtools view -bS -o {align_dir}/${{nm}}.bam {align_dir}/${{nm}}.sam\n'
                       'rm {align_dir}/${{nm}}.sam\n')

	align_template = ('#!/bin/bash -l\n'
    		           '#SBATCH -A b2012025\n'
    		           '#SBATCH -J {sam}_align\n'
    		           '#SBATCH -p core -n 2 \n'	
    		           '#SBATCH -t	10:00:00\n'
    		           '#SBATCH --mail-type=FAIL\n'
    		           '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
    		           'module load bioinfo-tools\n'
    		           'module load samtools/1.3\n'
			   ''+align_module+''
		 	   'cd '+proj_dir+'\n'
			   'if [[ $(ls {sam}/Rawdata/*gz | wc -l) -gt 0 ]]; then gzip -d {sam}/Rawdata/*gz; fi\n'
			   'for run in {fq_files};do\n'
			   'nm=$(basename ${{run}})\n'
			   'nm=${{nm/.fastq/}}\n'
			   ''+align_block+''
			   'done\n'
                        )
	samples = find_samples(proj_dir)
	for sam in samples.keys():
        	fq_fls = samples[sam]
        	sam_dir = os.path.join(proj_dir, sam)
		src_dir = os.path.join(sam_dir, 'scripts')
        	if not os.path.exists(src_dir):
        		os.mkdir(src_dir)
		align_dir = os.path.join(sam_dir,aligner)
		if not os.path.exists(align_dir):
			os.mkdir(align_dir)
        	job_file = os.path.join(src_dir, "{}_{}.sh".format(sam,aligner))
		with open(job_file, 'w') as jb_fl:
        		jb_fl.write(align_template.format(sam=sam, sam_dir=sam_dir,align_dir=align_dir,proj_dir=proj_dir,align_index=align_index, fq_files=" ".join(fq_fls)))
		#subprocess.check_call(['sbatch',job_file])
