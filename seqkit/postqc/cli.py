#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.postqc import postqc as pq

@click.group()
def postqc():
    """ postqc methods entry point """
    pass

## postqc subcommands
@postqc.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform postqc')
@click.option('--genefile', required=True, type=click.STRING,help='UCSC gene region file')
@click.option('-s','--sample',type=click.STRING,help='Sample to run')


@click.pass_context
def postqc(ctx, project, genefile, sample):
    """ Commands to run deepTools QC """
    pq.bamcov(project, genefile, sample)

