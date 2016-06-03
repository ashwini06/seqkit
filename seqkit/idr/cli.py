#!/usr/bin/env python
""" CLI for the idr subcommand
"""
import click
from seqkit.idr import idr as ir

@click.group()
def idr():
        """ IDR_analysis methods entry point """
        pass

## Peakcalling subcommands
@idr.command()
@click.option('-p','--project', required=True, type=click.STRING, help='Project to perform alignment')
@click.option('-i','--input_file', required=True, type=click.STRING, help='Input file containing which files to perform idr analysis between replicates')
@click.pass_context
def idr(ctx,project,input_file):
    """ IDR_analysis methods to check consistence across biological replicates"""
    ir.run_idr(project,input_file)










