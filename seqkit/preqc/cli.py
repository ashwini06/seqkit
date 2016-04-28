#!/usr/bin/env python
""" CLI for the QC subcommand
"""
import click
from seqkit.preqc import preqc as qc
        
@click.command()
@click.option('-p', '--project', type=click.STRING, help='Project to perform preQC')
@click.option('-o' , '--out_Dir', type=click.STRING, help='Path to output folder')
@click.pass_context
def preqc(ctx, project):
	""" pre QC methods and utilities """
	qc.run_qc("/proj/b2012025/RAW_DATA/ChIP_histone/lizzy_andersson_ChIP_ESC_wt_57-64")
	cmd = ['module','load','FastQC/0.11.5']
	subprocess.check_all(cmd)
	cmd = ['fastqc', samples, '-o' ]
	subprocess.check_all(cmd)
