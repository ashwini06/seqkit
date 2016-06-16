#!/usr/bin/env python
""" CLI for the peakanalysis subcommand
"""
import click
from seqkit.peakcall import peakanalysis as pa


@click.group()
def peakanalysis():
    """ Peakanalysis method entry point """
    pass

## peakcall subcommand
@peakanalysis.command()
@click.option('-p','--project', required=True, type=click.STRING, help='Project to perform peak-calling')
@click.option('-i','--input_file', required=True, type=click.STRING, help='Input file containing which files to be used for treatmentVsControl')
@click.option('-m','--mode', default='TF', type=click.STRING, help='type of ChIP-data (Transcription factor(TF) or Histone modifications(HM), default is TF)')
@click.option('--peak_call', default='macs2', type=click.STRING, help='which peak-calling to use macs2/sissrs/danpos2, default is "macs2"')
@click.option('--peakannotate', is_flag=True, help="do peak annoanotation")
def peakcall(project,input_file,mode,peak_call,peakannotate):
	""" Peak calling methods to call the ChIPSeq signals"""
	pa.run_peakcall(project,input_file,mode,peak_call,peakannotate)

@peakanalysis.command()
@click.option('-p','--project', required=True, type=click.STRING, help='Project to perform preak-annotation')
@click.option('--peak_call', default='macs2', type=click.STRING, help='which peak-calling to use macs2/sissrs/danpos2, default is "macs2"')
@click.option('--slurm', is_flag=True, help='Run the conversion as slurm jobs')
@click.pass_context
def peakanno(ctx,project,peak_call,slurm):
    """ Commands to annotate peak-called regions """
    pa.run_peakanno(project,peak_call,slurm)

