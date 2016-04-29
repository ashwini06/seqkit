#!/usr/bin/env python
""" CLI for the QC subcommand
"""
import click
import subprocess
from seqkit.preqc import preqc as qc
        
@click.command()
@click.option('-p', '--project', required=True, type=click.STRING, help='Project to perform preQC')
@click.option('-o' , '--out_dir', type=click.STRING, help='Path to output folder')
@click.pass_context
def preqc(ctx, project, out_dir):
	""" pre QC methods and utilities """
	qc.run_qc(project)

