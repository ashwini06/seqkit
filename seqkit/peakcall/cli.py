#!/usr/bin/env python
""" CLI for the peakanalysis subcommand
"""
import click
from seqkit.peakanalysis import peakanalysis as pa

@click.group()
def peakanalysis():
    """ Peakanalysis methods entry point """
    pass

## peakcall subcommand
@peakanalysis.command()
@click.option('-p','--project', required=True, type=click.STRING, help='Project to perform peak-calling')
@click.option('-i','--input_file', required=True, type=click.STRING, help='Input file containing which files to be used for treatmentVsControl')
@click.option('-m','--mode', default='TF', type=click.STRING, help='type of ChIP-data (Transcription factor(TF) or Histone modifications(HM), default is TF)')
@click.option('--peak_call', default='macs2', type=click.STRING, help='which peak-calling to use macs2/sissrs/danpos2, default is "macs2"')
def peakcall(project,input_file,mode,peak_call):
	""" Peak calling methods to call the ChIPSeq signals"""
	pa.run_peakcall(project,input_file,mode,peak_call)
