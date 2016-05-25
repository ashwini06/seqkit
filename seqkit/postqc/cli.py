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
@click.option('-i','--input_file',required=True,type=click.STRING,help='Input file containing which files to be used for treatment vs control')

@click.pass_context
def postqc(ctx, project, genefile, input_file):
    """ Commands to run deepTools QC """
    pq.bamcov(project, genefile, input_file)

