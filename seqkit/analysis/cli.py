#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.analysis import analysis as als
        
@click.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform bowtie alignment')
@click.option('-o', '--out_dir', type=click.STRING, help='Path to save outputs')
@click.pass_context
def analysis(ctx, project, out_dir):
	""" Commands to run analysis/pipelines """
	als.run_bowtie(project)
