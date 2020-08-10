"""
cdn Command Line Interface
author: Kevin Menden, DZNE Tübingen
This is the main file for executing the cdn program.
"""

# imports
import click
from run_deconccn import deconvolution

@click.group()
@click.version_option('1.0.0')
def cli():
    pass



"""
Processing mode
"""
@cli.command()
@click.argument(
    'ref_path',
    type = click.Path(exists=True),
    required = True
)
@click.argument(
    'mix_path',
    type = click.Path(exists=True),
    required = True
)
@click.argument(
    'marker_path',
    type = click.Path(exists=True),
    required = True
)

@click.option(
    '--unknown',
    default = False,
    help = 'unknown content'
)
def deconccnEstimate(ref_path, mix_path, marker_path, unknown):
    """ DeconCCN run """
    deconvolution(ref_path, mix_path, marker_path, unknown)

def main():
    text = """
     ___    ____    __    ___   _    _   __    __  _    _
    |   \  |       /     /   \  |\   |  /     /    |\   |   
    |    \ |____  /     |     | | \  | /     /     | \  | 
    |    / |      \     |     | |  \ | \     \     |  \ | 
    |___/  |____   \__   \___/_ |   \|  \__   \__  |   \| 
    """
    click.echo(click.style(text, fg='red'))
    cli()