#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.analysis import analysis as als
        
@click.command()
@click.option('-o', '--output', type=click.STRING, help='Path to save outputs')
@click.pass_context
def analysis(ctx, output):
	""" Commands to run analysis/pipelines """
	als.run_analysis(output)
