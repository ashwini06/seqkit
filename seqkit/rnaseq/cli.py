#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.rnaseq import rnaseq as rs
        
@click.group()
def rnaseq():
    """ rnaseq methods entry point """
    pass

## analysis subcommands
@rnaseq.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform alignment')
@click.option('-a','--aligner', default='STAR',type=click.STRING,help='which aligner to use STAR/Tophat, default is "STAR"')
@click.option('-s','--sample',type=click.STRING,help='Sample to run')
@click.pass_context
def htcuff(ctx, project, aligner, sample):
    """ Commands to run  htseq-Count and Cufflinks """ 
    rs.run_htcuff(project, aligner, sample)
