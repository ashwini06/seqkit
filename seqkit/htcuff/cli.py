#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.htcuff import htcuff as hc
        
## htcuff command 
@click.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform alignment')
@click.option('-a','--aligner', default='STAR',type=click.STRING,help='which aligner to use STAR/Tophat, default is "STAR"')
@click.option('-s','--sample',type=click.STRING,help='Sample to run')
@click.pass_context
def htcuff(ctx, project, aligner, sample):
    """ Commands to run  htseq-Count and Cufflinks """ 
    hc.run_htcuff(project, aligner, sample)
