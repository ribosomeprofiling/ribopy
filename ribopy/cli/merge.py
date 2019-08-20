# -*- coding: utf-8 -*-

from .main import *
from ..merge import merge_ribo_files

@cli.command()
@click.argument(
    "out_ribo_path",
    nargs = 1,
    type  = click.Path(exists = False))
@click.argument(
    "in_ribo_paths",
    nargs = -1,
    type  = click.Path(exists = True))
def merge(out_ribo_path, in_ribo_paths):
    """
    Merges a set of given ribo files into one ribo file.
                     
    The input ribo files are merged into a new ribo file.
    The resulting ribo file has the union of the experiments
    of the input files

    The ribo files to be merged must be compatible:
    They must have the same\n 
    1) Reference Name \n
    2) Transcript Names & Transcript Lengths \n
    3) Annotation \n
    4) Ribo file parameters: \n
       a) Metagene Radius  \n
       b) Left Span & Right Span \n
       c) Min & Max Read Length \n
    5) They can not have overlapping experiment names

    This option is especially useful for Ribo-Seq Data
    processed together and needs to be analyzed together.
    """
    
    merge_ribo_files(out_ribo_path, in_ribo_paths)
