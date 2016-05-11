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
@click.option('-p','--project',required=True, type=click.STRING,help='Project to perform alignment')
@click.option('-o', '--out_dir', type=click.STRING, help='Path to save outputs')
@click.option('-pc','--peak_call', default='macs2',type=click.STRING,help='which peak-calling to use macs2/sissrs/danpos2, default is "macs2"')
@click.pass_context
def peakcall(ctx,project,out_dir,peak_call):
	""" Peak calling methods to call the ChIPSeq signals"""
	pc.run_peakcall(project,peak_call)
