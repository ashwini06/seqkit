#!/usr/bin/env python
""" CLI for the motifsearch subcommand
"""
import click
from seqkit.motif import motif as ms


@click.group()
def motif():
    """ Motif analysis method entry point """
    pass


## motifanalysis subcommand
@motif.command()
@click.option('-p','--project', required=True, type=click.STRING, help='Project to perform motif-analysis')
@click.option('--peak_call', default='macs2', type=click.STRING, help='peak-calling folder to use macs2/sissrs/danpos2, default is "macs2"')
@click.pass_context
def denovo(ctx,project,peak_call):
    """ De-novo motif analysis """
    ms.run_denovo(project,peak_call)


