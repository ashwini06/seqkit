#!/usr/bin/env python
""" CLI for the QC subcommand
"""
import click
from seqkit.preqc import preqc as qc
        
@click.command()
@click.option('-o', '--output', type=click.STRING, help='Path to save outputs')
@click.pass_context
def preqc(ctx, output):
	""" pre QC methods and utilities """
	qc.run_qc(output)
