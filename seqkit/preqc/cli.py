#!/usr/bin/env python
""" CLI for the QC subcommand
"""
import click
from seqkit.preqc import preqc as qc
        
@click.command()
@click.option('-p', '--project', type=click.STRING, help='Project to perform preQC')
@click.pass_context
def preqc(ctx, project):
	""" pre QC methods and utilities """
	qc.run_qc(output)
