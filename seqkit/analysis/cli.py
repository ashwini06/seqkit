#!/usr/bin/env python
""" CLI for the analysis subcommand
"""
import click
from seqkit.analysis import analysis as als
        
@click.group()
def analysis():
	""" Analysis methods entry point """
	pass

## analysis subcommands
@analysis.command()
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform alignment')
@click.option('-a','--aligner', default='bowtie2',type=click.STRING,help='which aligner to use bowtie2/bwa, default is "bowtie2"')
@click.option('--bam_to_bed', is_flag=True, help="Convert the aligned bam files to bed files")
@click.pass_context
def align(ctx, project, aligner, bam_to_bed):
	""" Commands to run analysis/pipelines """
	als.run_align(project, aligner, bam_to_bed)

@analysis.command()
@click.option('-p','--project',required=True, type=click.STRING,help='project to perform bam2bed conversion')
@click.option('-a','--aligner', default='bowtie2',type=click.STRING,help='which aligner was used to align, default is "bowtie2"')
@click.option('-s', '--slurm', is_flag=True, help='Run the conversion as slurm jobs')
@click.pass_context
def bamTobed(ctx, project, aligner, slurm):
	""" Commands to convert bam files to bed files """
	als.run_b2b(project, aligner, slrum)
