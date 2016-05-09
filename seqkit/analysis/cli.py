#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.analysis import analysis as als
        
@click.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform bowtie alignment')
@click.option('-o', '--out_dir', type=click.STRING, help='Path to save outputs')
@click.option('-a','--aligner', default='bowtie2',type=click.STRING,help='which aligner to use bowtie2/bwa, default is "bowtie2"')
@click.pass_context
def analysis(ctx, project, out_dir, aligner):
	""" Commands to run analysis/pipelines """
	als.run_align(project,aligner)
