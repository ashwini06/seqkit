#!/usr/bin/env python
"""Methods and functionalites to do analysis"""
import os
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def run_b2b(project, aligner, sample=None, slurm=False, job_file=None):
    """ Will run the bam to bed file conversion """
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)

    if sample:
        if os.path.isdir(os.path.join(proj_dir, sample)):
            samples = [sample]
        else:
            raise SystemExit("Given sample {} is not found in project directory {}".format(sample, proj_dir))
    else:
        samples = find_samples(proj_dir)

    for sam in samples:
        sbatch_template = ('#!/bin/bash -l\n'
                   '#SBATCH -A b2012025\n'
                   '#SBATCH -J {sam}_bam2bed\n'
                   '#SBATCH -p core -n 1 \n'    
                   '#SBATCH -t 3:00:00\n'
                   '#SBATCH --mail-type=FAIL\n'
                   '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n')
        template_b2b = ('## run bam to bed\n'
                'module load BEDTools/2.11.2\n'
                'for bam in $(ls --color=never {sam_dir}/alignment_{aligner}/bam_files/*sorted_rmdup.bam);do\n'
                'bed_fl=${{bam/.bam/.bed}}\n'
                'bed_fl=${{bed_fl/bam_files/bedfiles}}\n'
                'bed_uniq_fl=${{bed_fl/.bed/_uniq.bed}}\n'
                'bamToBed -i ${{bam}} > ${{bed_fl}}\n'
                'awk -F\\\\t -v \'OFS=\\t\' \'{{print "chr"$1,$2,$3,".",$5,$6}}\' ${{bed_fl}} | sort -u > ${{bed_uniq_fl}}\n'
                'rm ${{bed_fl}}\n'
                'done\n')
        sam_dir = os.path.join(proj_dir, sam)
        bed_dir = os.path.join(sam_dir, "alignment_{}".format(aligner), "bedfiles")
        if not os.path.exists(bed_dir): 
            os.mkdir(bed_dir)
        if job_file:
            with open(job_file, 'a') as jb_fl:
                jb_fl.write(template_b2b.format(sam_dir=sam_dir, aligner=aligner))
            return
        if slurm:
            job_file = os.path.join(sam_dir,"scripts","{}_{}_bamTobed.sh".format(sam, aligner))
            template_b2b = sbatch_template + template_b2b
            with open(job_file, 'w') as jb_fl:
                jb_fl.write(template_b2b.format(sam=sam, sam_dir=sam_dir, aligner=aligner))
            subprocess.check_call(['sbatch',job_file])
            job_file = None


def run_align(project, aligner, sample, bam_to_bed):
    """Will run the preferred-alignment"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    bed_dir = ''
    if aligner == "bwa":
        align_module = 'module load bwa/0.7.12\n'
        align_index = '/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa'
        align_block = ('bwa aln {align_index} ${{fq}} > {align_dir}/${{nam}}.sai\n'
                       'bwa samse {align_index} {align_dir}/${{nam}}.sai ${{fq}} | samtools view -Sb - > {align_dir}/${{nam}}.bam\n'
                       'rm ${{nam}}.sai\n')
    elif aligner == "bowtie2":
        align_module = 'module load bowtie2/2.2.6\n'
        align_index = "/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome"
        align_block =  ('bowtie2 -t -p 8 -k2 --very-sensitive -x {align_index} -q ${{fq}} -S {align_dir}/${{nam}}.sam 2> {align_dir}/${{nam}}_bowtie2.log\n\n'
                        'samtools view -bS -o {align_dir}/${{nam}}.bam {align_dir}/${{nam}}.sam\n\n'
                        'rm {align_dir}/${{nam}}.sam\n\n')
    elif aligner == "bowtie":
        align_module = 'module load bowtie/1.1.2\n'
        align_index = "/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Sequence/BowtieIndex/genome"
        align_block =  ('bowtie -q -m 1 -v 3 --best --strata {align_index} ${{fq}} -S {align_dir}/${{nam}}.sam 2>{align_dir}/${{nam}}_bowtie.log\n\n'
                        'samtools view -bS -o {align_dir}/${{nam}}.bam {align_dir}/${{nam}}.sam\n\n')

    
    elif aligner == "STAR":
        align_module = 'module load star/2.3.1o\n'
        align_index = "/pica/data/uppnex/igenomes_new/Mus_musculus/Ensembl/GRCm38/Sequence/STARIndex"
        align_block = ("STAR --genomeDir {align_index} --readFilesIn ${{fq}} --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {align_dir}/${{nam}} --outSAMmode Full --runThreadN 8 --outFilterType BySJout --alignSJDBoverhangMin 1 --outFilterMismatchNmax 5\n\n")
        
    align_template = ('#!/bin/bash -l\n'
                      '#SBATCH -A b2012025\n'
                      '#SBATCH -J {sam}_align\n'
                      '#SBATCH -p core -n 4 \n' 
                      '#SBATCH -t 10:00:00\n'
                      '#SBATCH --mail-type=FAIL\n'
                      '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                      '#SBATCH -e {sam_dir}/scripts/{sam}_align.stderr\n'
                      '#SBATCH -o {sam_dir}/scripts/{sam}_align.stdout\n'
                      'module load bioinfo-tools\n'
                      'module load samtools/1.3\n'
                      ''+align_module+''
                      'if [[ $(ls --color=never {sam_dir}/Rawdata/*.gz | wc -l) -gt 0 ]]; then gzip -d {sam_dir}/Rawdata/*.gz; fi\n'
                      'if [[ $(ls --color=never {sam_dir}/Rawdata/*zip | wc -l) -gt 0 ]]; then unzip {sam_dir}/Rawdata/*zip; fi\n'
                      'if [[ $(ls --color=never {sam_dir}/Rawdata/*gz | wc -l) -gt 0 ]]; then gzip -d {sam_dir}/Rawdata/*gz; fi\n\n'
                      'for fq in $(ls --color=never {sam_dir}/Rawdata/*.fastq);do\n'
                      'nm=$(basename ${{fq}})\n'
                      'nm=${{nm/_*/}}\n' 
                      'nam="{sam}_"${{nm}}\n\n'
                      ''+align_block+''
                      'samtools view -H {align_dir}/${{nam}}.bam | sed -e \'s/SN:\([0-9XY]\)/SN:chr\\1/\' -e \'s/SN:M/SN:chrM/\' | samtools reheader - {align_dir}/${{nam}}.bam > {align_dir}/${{nam}}_v1.bam\n\n'
                      'mv {align_dir}/${{nam}}_v1.bam {align_dir}/${{nam}}.bam\n\n'
                      'samtools sort -T {align_dir}/temp -o {align_dir}/${{nam}}_sorted.bam {align_dir}/${{nam}}.bam\n\n'
                      'java -jar /pica/sw/apps/bioinfo/picard/1.92/milou/MarkDuplicates.jar INPUT={align_dir}/${{nam}}_sorted.bam OUTPUT={align_dir}/${{nam}}_sorted_rmdup.bam METRICS_FILE={align_dir}/${{nam}}_picardmetrics.txt REMOVE_DUPLICATES=True\n\n'
                      'samtools index {align_dir}/${{nam}}_sorted_rmdup.bam\n\n'
                      'samtools index {align_dir}/${{nam}}_sorted.bam\n\n'
                      'rm {align_dir}/${{nam}}.bam\n\n'
                      'rm {align_dir}/${{nam}}.sam\n\n'
                      'done\n')

    if sample:
        if os.path.isdir(os.path.join(proj_dir, sample)):
            samples = [sample]
        else:
            raise SystemExit("Given sample {} is not found in project directory {}".format(sample, proj_dir))
    else:
        samples = find_samples(proj_dir)

    for sam in samples:
        sam_dir = os.path.join(proj_dir, sam)
        src_dir = os.path.join(sam_dir, 'scripts')
        if not os.path.exists(src_dir):
            os.mkdir(src_dir)
        align_dir = os.path.join(sam_dir,"alignment_{}".format(aligner),"bam_files")
        if not os.path.exists(align_dir):
            os.makedirs(align_dir)
        job_file = os.path.join(src_dir, "{}_{}.sh".format(sam,aligner))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(align_template.format(sam=sam, sam_dir=sam_dir, align_dir=align_dir, align_index=align_index))
        if bam_to_bed:
            run_b2b(project=project, aligner=aligner, slurm=True, sample=sam, job_file=job_file)
        subprocess.check_call(['sbatch',job_file])
