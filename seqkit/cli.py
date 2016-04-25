#!/usr/bin/env python
import os
import logging

from pkg_resources import iter_entry_points

import click

from seqkit import __version__
from seqkit import load_yaml_config

logger = logging.getLogger(__name__)


@click.group()
@click.version_option(__version__)
# Priority for the configuration file is: environment variable > -c option > default
@click.option('-c', '--config-file',
			  default=os.path.join(os.environ['HOME'], '.seqkit/seqkit.yaml'),
			  envvar='SEQKIT_CONFIG',
			  type=click.File('r'),
			  help='Path to SeqKit configuration file')
@click.pass_context
def main(ctx, config_file):
	""" Toolkit for the Automation of QC and Analyses """
	ctx.obj = {}
	config = load_yaml_config(config_file)
	logger.debug('starting up CLI')

#Add subcommands dynamically to the CLI
for entry_point in iter_entry_points('seqkit.subcommands'):
	main.add_command(entry_point.load())