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
@click.option('-a','--aligner', default='bowtie',type=click.STRING,help='which aligner to use bowtie/bowtie2/bwa, default is "bowtie"')
@click.option('-s','--sample',type=click.STRING,help='Sample to run')
@click.option('--bam_to_bed', is_flag=True, help="Convert the aligned bam files to bed files")
@click.pass_context
def align(ctx, project, aligner, sample, bam_to_bed):
	""" Commands to run analysis/pipelines """
	als.run_align(project, aligner, sample, bam_to_bed)

@analysis.command()
@click.option('-p','--project',required=True, type=click.STRING,help='project to perform bam2bed conversion')
@click.option('-a','--aligner', default='bowtie',type=click.STRING,help='which aligner was used to align, default is "bowtie"')
@click.option('--slurm', is_flag=True, help='Run the conversion as slurm jobs')
@click.option('-s','--sample',type=click.STRING,help='Sample to run')
@click.pass_context
def bamTobed(ctx, project, aligner, sample, slurm):
	""" Commands to convert bam files to bed files """
	als.run_b2b(project, aligner, sample, slurm)
