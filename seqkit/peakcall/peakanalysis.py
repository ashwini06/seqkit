#!/usr/bin/env python
""" Methods and functionalites to perform peak-calling """
import os
import subprocess
from glob import glob
from seqkit import CONFIG as conf
from seqkit.utils import col_match 


def run_peakanno(project,peak_call,slurm=False,job_file=None):
    """ Will run the peak annotation on the peak called regions """
    root_dir = conf.get('root_dir','')
    proj_dir = os.path.join(root_dir,project)
    samples = map(str,glob(os.path.join(proj_dir, "*", "*{}*".format(peak_call))))
    TSS_cmd = conf.get ('anno_TSS','')
    NDG_cmd = conf.get ('anno_NDG','') 
    sbatch_template = ('#!/bin/bash -l\n'
    '#SBATCH -A b2012025\n'
    '#SBATCH -J {nm}_peakanno\n'
    '#SBATCH -p core -n 1 \n'
    '#SBATCH -t 3:00:00\n'
    '#SBATCH --mail-type=FAIL\n'
    '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n')
    template_peakanno = ('## Running peak-annotations\n'
    'for bed in $(ls --color=never {peaks_dir}/*narrowPeak);do\n'
    'cut -f1-6 $bed > {peaks_dir}/{nm}_annotate \n'
    ''+TSS_cmd+'\n'
    ''+NDG_cmd+'\n'
    'python '+col_match.__file__+'{peaks_dir}/{nm}_annotate.tss {peaks_dir}/{nm}_annotate.ndg comb_fl "merge"\n'
    'rm {peaks_dir}/{nm}_annotate\n'
    'done\n')             
    
#    colmatch({peaks_dir}/{nm}_annotate.tss,{peaks_dir}/{nm}_annotate.ndg,comb_fl,"merge")
 
    for sam in samples:
        sam_dir=os.path.split(sam)[0]
        nm=os.path.basename(sam_dir)
        annotate_dir = os.path.join(sam_dir,"peakannotate")
        if not os.path.exists(annotate_dir):
            os.makedirs(annotate_dir)
        if job_file:
            with open(job_file,'a') as jb_fl:
                jb_fl.write(template_peakanno.format(peaks_dir=sam,nm=nm,anno_dir=annotate_dir))
            return
        if slurm:
            job_file = os.path.join(sam_dir,"scripts","{}_peakannotate.sh".format(nm))
            template_anno = sbatch_template+template_peakanno
            with open(job_file,'w') as jb_fl:
                jb_fl.write(template_anno.format(peaks_dir=sam,nm=nm,anno_dir=annotate_dir))
            job_file=None



def run_peakcall(project, input_file, mode, peak_call,peakannotate):
    """ Will run the prefered peak-calling software """
    root_dir = conf.get('root_dir','')  
    proj_dir = os.path.join (root_dir,project)
    load_module = ('module load MACS/2.1.0\n')
    sbatch_template = ('#!/bin/bash -l\n'
                       '#SBATCH -A b2012025\n'
                       '#SBATCH -J {name}_peakcall\n'
                       '#SBATCH -p core -n 1 \n'
                       '#SBATCH -t 5:00:00\n'
                       '#SBATCH -o '+proj_dir+'/{treat}/scripts/{name}_peakcall.stdout\n'
                       '#SBATCH -e '+proj_dir+'/{treat}/scripts/{name}_peakcall.stderr\n'
                       '#SBATCH --mail-type=FAIL\n'
                       '#SBATCH --mail-user=\'ashwini.jeggari@scilifelab.se\'\n\n'
                       'module load bioinfo-tools\n'
                        )
    if mode == "TF":
        if peak_call == "macs2":
            macs2_cmd = conf.get('macs2_TF','')
            template = ('## Running Peak-calling for TF-ChIP data\n'
                        ''+load_module+''
            ''+macs2_cmd+''
                )
        else:
            raise SystemExit("Please mention the type of peak caller - macs2")
    elif mode == "HM":   
        if peak_call == "macs2":
            macs2_cmd = conf.get('macs2_HM','')
            template = ('# Running macs2 peak-calling for HM data\n'
                        ''+load_module+''
                        ''+macs2_cmd+'')
    
        elif peak_call == "danpos2":
            danpos_path = "cd /home/ashwini/softwares/danpos-2.2.2"
            danpos_cmd = conf.get('danpos2_dpeak','')
            template = ('# Running danpos2 peakcalling for HM data\n'
                        ''+danpos_cmd+'')
        else:
            raise SystemExit("Please mention the type of peak_Caller (macs2/danpos2)")
    else:
        raise SystemExit("Please mention the type of mode - either TF or HM")

        
    pk_file = open(input_file,'r')
    pk_file.next()
    for ln in iter(pk_file):    
        ln = ln.strip()
        ln =  ln.split('\t')
        treat = ln[0]
        ctrl = ln[1]
    treat_fl = glob("{}/{}/alignment_*/bedfiles/{}*rmdup_uniq.bed".format(proj_dir,treat,treat))
    control_fl = glob("{}/{}/alignment_*/bedfiles/{}*rmdup_uniq.bed".format(proj_dir,ctrl,ctrl))
        peaks_dir = os.path.join(proj_dir,treat,"{}_{}".format(peak_call,mode))
        if not os.path.exists(peaks_dir):
            os.makedirs(peaks_dir)
        for sam in treat_fl:
            suf_s = os.path.basename(sam)
            suf_s = suf_s.replace("_sorted_rmdup_uniq.bed","")
            for con in control_fl:
                con_c = os.path.basename(con)
                con_c = con_c.replace("_sorted_rmdup_uniq.bed","")
                name = "{}_Vs_{}".format(suf_s,con_c)
                job_fl = os.path.join(proj_dir,treat,"scripts","{}_peakcall.sh".format(name))
                template_pc = sbatch_template + template
                with open(job_fl,'w') as jb_fl:
                    jb_fl.write(template_pc.format(name=name,treat=treat, treatment=sam, control=con,peaks_dir=peaks_dir))
                if peakannotate:
                    run_peakanno(project=project,peak_call=peak_call,slurm=True,job_file=job_fl)        
                #subprocess.check_call(['sbatch',job_fl]) 

