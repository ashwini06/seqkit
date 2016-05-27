#!/usr/bin/env python
""" CLI for the peakcall subcommand
"""
import click
from seqkit.peakcall import peakcall as pc

@click.group()
def peakcall():
        """ Peakcalling methods entry point """
        pass

## Peakcalling subcommands
@peakcall.command()
@click.option('-p','--project', required=True, type=click.STRING, help='Project to perform alignment')
@click.option('-i','--input_file', required=True, type=click.STRING, help='Input file containing which files to be used for treatmentVsControl')
@click.option('-m','--mode', default='TF', type=click.STRING, help='type of ChIP-data (Transcription factor(TF) or Histone modifications(HM), default is TF)')
@click.option('--peak_call', default='macs2', type=click.STRING, help='which peak-calling to use macs2/sissrs/danpos2, default is "macs2"')
@click.pass_context
def peakcall(ctx,project,input_file,mode,peak_call):
	""" Peak calling methods to call the ChIPSeq signals"""
	pc.run_peakcall(project,input_file,mode,peak_call)
