#!/usr/bin/env python
"""Methods and functionalites to do analysis"""
import os
import pdb
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def run_b2b(project, aligner, slurm=False, samples=None, job_file=None):
	""" Will run the bam to bed file conversion """
	root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
	template_b2b = ('## run bam to bed\n'
		    'module load BEDTools/2.11.2\n'
		    'bamToBed -i {bam_file} > {bed_file}\n'
	        'awk -F\\\\t -v \'OFS=\\t\' \'{{print chr$1,$2,$3,".",$5,$6}}\' {bed_file} | sort -u > {bed_uniq_file}\n'
	 	    'rm {bed_file}\n')
    if not samples:
        samples = find_samples(proj_dir, file_type="bam", aligner=aligner):

    for sam in samples:
        bam_fls = samples[sam]
		sam_dir = os.path.join(proj_dir, sam)
		bed_dir = os.path.join(sam_dir, "alignment_{}".format(aligner), "bedfiles")
		if not os.path.exists(bed_dir):	
			os.mkdir(bed_dir)
        if slurm:
            if not job_file:
		        job_file = os.path.join(sam_dir,"scripts","{}_{}_bamTobed.sh".format(sam, aligner))
		    with open(job_file, 'a') as jb_fl:
	       	    jb_fl.write(template_b2b.format(bam_file=bam,bed_file=bed_file,bed_uniq_file=bed_uniq_file))

def run_align(project, aligner, bam_to_bed):
	"""Will run the preferred-alignment"""
	root_dir = conf.get('root_dir','')
	proj_dir = os.path.join(root_dir, project)
    bed_dir = ''
	if aligner == "bwa":
		align_module = 'module load bwa/0.7.12\n'
		align_index = '/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa'
    	align_block = ('bwa aln {align_index} ${{run}} > {align_dir}/${{nam}}.sai\n'
                       'bwa samse {align_index} {align_dir}/${{nam}}.sai ${{run}}.fastq | samtools view -Sb - > {align_dir}/${{nam}}.bam\n'
				       'rm ${{nam}}.sai\n')
	elif aligner == "bowtie2":
		align_module = 'module load bowtie2/2.2.6\n'
		align_index = "/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome"
		align_block =  ('bowtie2 -t -p 8 -k2 --very-sensitive -x {align_index} -q ${{run}} -S {align_dir}/${{nam}}.sam > {align_dir}/${{nam}}_bowtie2.log\n'
                        'samtools view -bS -o {align_dir}/${{nam}}.bam {align_dir}/${{nam}}.sam\n'
		                'rm {align_dir}/${{nam}}.sam\n')
	
	align_template = ('#!/bin/bash -l\n'
    		           '#SBATCH -A b2012025\n'
    		           '#SBATCH -J {sam}_align\n'
    		           '#SBATCH -p core -n 2 \n'	
    		           '#SBATCH -t	10:00:00\n'
    		           '#SBATCH --mail-type=FAIL\n'
    		           '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
    		           'module load bioinfo-tools\n'
    		           'module load samtools/1.3\n'
			           +align_module+
		 	           'cd '+proj_dir+'\n'
			           'if [[ $(ls {sam}/Rawdata/*gz | wc -l) -gt 0 ]]; then gzip -d {sam}/Rawdata/*gz; fi\n'
			           'for run in {fq_files};do\n'
			           'nm=$(basename ${{run}})\n'
#			           'nm=${{nm/.fastq/}}\n'
			           'nm=${{nm/_*/}}\n' 
			           'nam={sam}"_"${{nm}}\n'
			           +align_block+
			           'samtools sort -m 10000000 {align_dir}/${{nam}}.bam {align_dir}/${{nam}}_sorted'
                       'samtools index {align_dir}/${{nam}}_sorted'
			           'rm {align_dir}/${{nam}}.bam\n'
			           'done\n')

	samples = find_samples(proj_dir)
	for sam in samples.keys():
        fq_fls = samples[sam]
        sam_dir = os.path.join(proj_dir, sam)
		src_dir = os.path.join(sam_dir, 'scripts')
        if not os.path.exists(src_dir):
        	os.mkdir(src_dir)
		align_dir = os.path.join(sam_dir,"alignment_{}".format(aligner),"bam_files")
		if not os.path.exists(align_dir):
			os.makedirs(align_dir)
        job_file = os.path.join(src_dir, "{}_{}.sh".format(sam,aligner))
		with open(job_file, 'w') as jb_fl:
        	jb_fl.write(align_template.format(sam=sam, sam_dir=sam_dir,align_dir=align_dir,proj_dir=proj_dir,
                                              align_index=align_index, fq_files=" ".join(fq_fls), bed_dir=bed_dir))
        if bam_to_bed:
            tmp_bam_file_names = []
            for fq in fq_fls:
                f_nm = os.path.join(align_dir, "{}_sorted.bam".format(os.path.basename(fq).replace('.gz','').replace('.fastq','')))
                tmp_bam_file_names.append(f_nm)
            run_b2b(project=project, aligner=aligner, slurm=True, samples={sam:tmp_bam_file_names}, job_file=job_file)
		subprocess.check_call(['sbatch',job_file])




	
