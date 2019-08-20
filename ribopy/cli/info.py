from .main import *
import click
from ..info import ribo_file_info
                   
###########################################################

@cli.command()
@click.argument('ribofile', type=click.Path(exists=True ))
def info(ribofile):
    """
    Displays a summary information about the given ribo file
    """

    display_str = ribo_file_info(ribo_file = ribofile)

    print( display_str)


