#!/usr/bin/env python
"""Methods and functionalites to do a pre QC"""
import os
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils.find_samples import find_samples

def bamcov(project, genefile, input_file, mode):
    """Will run the postqc"""
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir, project)
    if mode == "scale":
        assign_mode = conf.get('computematrix_scale','')
    else: 
        assign_mode = conf.get('computematrix_TSS','')
    sbatch_template = ('#!/bin/bash -l\n'
                       '#SBATCH -A b2012025\n'
                       '#SBATCH -J {name}_postqc\n'
                       '#SBATCH -p core -n 3 \n'
                       '#SBATCH -t 4:00:00\n'
                       '#SBATCH -e '+proj_dir+'/{sample}/scripts/{name}_postqc.stderr\n'
                       '#SBATCH -o '+proj_dir+'/{sample}/scripts/{name}_postqc.stdout\n'
                       '#SBATCH --mail-type=FAIL\n'
                       '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                       'module load bioinfo-tools\n'
                       'module load deepTools/2.2.3\n'
                       #'module load ngsplot/2.61\n\n'
                       )
                    
    template = ('bamCompare -b1 {treatment} -b2 {control} --binSize 25 --ratio log2 --scaleFactorsMethod "readCount" -o {postqc_dir}/{treat}_Vs_{ctrl}_log2ratio_readcount.bw --normalizeUsingRPKM\n'
                ''+assign_mode+'\n'
                'plotHeatmap -m {postqc_dir}/matrix.mat.gz -out {postqc_dir}/{treat}_Vs_{ctrl}_heatmap_v2.png --heatmapHeight 25 --heatmapWidth 3 --whatToShow \'heatmap and colorbar\' --sortUsing max\n'
                )


    bed_file = genefile
    pk_file = open(input_file,'r')
    pk_file.next()
    for ln in iter(pk_file):
        ln = ln.strip()
        ln =  ln.split('\t')
        treat = ln[0]
        ctrl = ln[1]
        postqc_dir = os.path.join(proj_dir,treat,"deepTools")
        if not os.path.exists(postqc_dir):
            os.mkdir(postqc_dir)
        treat_fl = glob("{}/{}/alignment_*/bam_files/{}*sorted_rmdup_v1.bam".format(proj_dir,treat,treat))
        control_fl = glob("{}/{}/alignment_*/bam_files/{}*sorted_rmdup_v1.bam".format(proj_dir,ctrl,ctrl))
        for sam in treat_fl:
            suf_s = os.path.basename(sam)
            suf_s = suf_s.replace("_sorted_rmdup_v1.bam","")
            for con in control_fl:  
                con_c = os.path.basename(con)
                con_c = con_c.replace("_sorted_rmdup_v1.bam","")
                name = "{}_Vs_{}".format(suf_s,con_c)
                job_file = os.path.join(proj_dir,treat,"{}/{}_{}.sh".format("scripts",name,"postqc"))
                template_pc = sbatch_template+template
                with open(job_file, 'w') as jb_fl:
                    jb_fl.write(template_pc.format(sample=treat, treat=suf_s, ctrl=con_c, name=name, treatment=sam, control=con, bed_file=bed_file, postqc_dir=postqc_dir))
#                subprocess.check_call(['sbatch',job_file])



                       #'bamCoverage --bam {treatment} --binSize 25 --normalizeUsingRPKM -o {postqc_dir}/{treat}_coverage.bw -bl /home/ashwini/mm10_blacklisted-regions.bed\n'
                       #'bamCoverage --bam {control} --binSize 25 --normalizeUsingRPKM -o {postqc_dir}/{ctrl}_coverage.bw -bl /home/ashwini/mm10_blacklisted-regions.bed\n'
                       #'plotFingerprint -b {treatment} {control} -plot {postqc_dir}/{treat}_Vs_{ctrl}_fingerprint.png --labels {treat} {ctrl}\n'
                       #'multiBamSummary bins --bamfiles '+proj_dir+'/*/alignment_bowtie/bam_files/*_rmdup.bam -out {postqc_dir}/results.npz\n'
                       #'plotCorrelation --corData {postqc_dir}/results.npz --plotFile {postqc_dir}/scatterplot.pdf --corMethod pearson --whatToPlot scatterplot --skipZeros\n'
                        #'ngs.plot.r -G mm10 -R genebody -C {treatment}:{control} -O {postqc_dir}/{treat}_Vs_{ctrl}.genebody -T {treat}\n'
                       #'ngs.plot.r -G mm10 -R tss -C {treatment}:{control} -O {postqc_dir}/{treat}_Vs_{ctrl}.tss -T {treat} -L 3000 -FL 3000\n'
                       #'if [ -e *cnt ]; then rm *cnt; fi\n'

