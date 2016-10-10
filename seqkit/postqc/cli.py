#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.postqc import postqc as pq

@click.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform postqc')
@click.option('--genefile', required=True, type=click.STRING,help='Bed file containing specified regions')
@click.option('-i','--input_file',required=True,type=click.STRING,help='Input file containing which files to be used for treatment vs control')
@click.option('-m','--mode',default='TSS',type=click.STRING,help='which method to run while calculating computeMatrix reference-point(TSS)/scale-regions(scale)')
@click.pass_context
def postqc(ctx,project, genefile, input_file, mode):
    """ Commands to run deepTools QC """
    pq.bamcov(project, genefile, input_file, mode)



#import click
#from seqkit.postqc import postqc as pq

#@click.group()
#def postqc():
#    """ postqc methods entry point """
#    pass

## postqc subcommands
#@postqc.command()
#@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform postqc')
#@click.option('--genefile', required=True, type=click.STRING,help='Bed file containing specified regions')
#@click.option('-i','--input_file',required=True,type=click.STRING,help='Input file containing which files to be used for treatment vs control')
#@click.option('-m','--mode',default='TSS',type=click.STRING,help='which method to run while calculating computeMatrix reference-point(TSS)/scale-regions(scale)')
#@click.pass_context
#def postqc(ctx,project, genefile, input_file, mode):
#    """ Commands to run deepTools QC """
#    pq.bamcov(project, genefile, input_file, mode)

