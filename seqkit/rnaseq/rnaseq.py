#!/usr/bin/env python
"""Methods and functionalites to do htseq and cufflinks"""
import os
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples


def run_htcuff(project, aligner, sample=None):
    """Will run the cuuflinks and htseq"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    align_template = ('#!/bin/bash -l\n'
                      '#SBATCH -A b2012025\n'
                      '#SBATCH -J {sam}_htcuff\n'
                      '#SBATCH -p core -n 10 \n' 
                      '#SBATCH -t 10:00:00\n'
                      '#SBATCH --mail-type=FAIL\n'
                      '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                      '#SBATCH -e {sam_dir}/scripts/{sam}_htcuff.stderr\n'
                      '#SBATCH -o {sam_dir}/scripts/{sam}_htcuff.stdout\n'
                      'module load bioinfo-tools\n'
                      'module load samtools/1.3\n'
                      'module load cufflinks/2.2.1\n'
                      'module load htseq/0.6.1\n\n'
                      'genome_fl=\"/pica/data/uppnex/igenomes/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf\"\n'
                      'for bam in $(ls --color=never {align_dir}/*_sorted*.bam);do\n'
                      'nm=$(basename ${{bam}})\n'
                      'nm=${{nm/.bam/}}\n'
                      'htseq-count -s reverse -q -f bam ${{bam}} ${{genome_fl}} > {ht_dir}/${{nm}}_counts.txt\n'
                      'cufflinks -p 8 --library-type fr-firststrand -G ${{genome_fl}} -o {cuff_dir}/${{nm}}_cufflinks ${{bam}}\n\n')


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
        align_dir = os.path.join(sam_dir,"alignment_{}".format(aligner),"bam_files")
        bam_fl =  ''.join(glob("{}/{}*"))
        ht_dir = os.path.join(sam_dir,'htseq')
        if not os.path.exists(ht_dir):
            os.mkdir(ht_dir) 
        cuff_dir = os.path.join(sam_dir,'cufflinks')
        if not os.path.exists(cuff_dir):
            os.mkdir(cuff_dir)           
        job_file = os.path.join(src_dir, "{}_{}.sh".format(sam,"htcuff"))
        with open(job_file, 'w') as jb_fl:
            jb_fl.write(align_template.format(sam=sam, sam_dir=sam_dir, align_dir=align_dir,ht_dir=ht_dir,cuff_dir=cuff_dir))
        #subprocess.check_call(['sbatch',job_file])
