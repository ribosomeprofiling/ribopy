# -*- coding: utf-8 -*-

import click

from .main import *

from ..plot import *
from ..settings import *



@cli.group()
def plot():
    """Generate some basic plots for ribo files."""
    pass


@plot.command()
@click.argument("ribo",      type=click.Path(exists=True))
@click.argument("experiments", nargs = -1,)
@click.option('-s', '--site',
              type     = click.Choice(['start', 'stop']),
              required = True,
              help     = 'Site type.')
@click.option('-o', '--out',
              required = True,
              help     = 'Output file.')
@click.option('--lowerlength',
              required = False ,
              type     = click.INT,
              help     = 'Lower read length')
@click.option('--upperlength',
              required = False ,
              type     = click.INT,
              help     ='Upper read length')
@click.option('-t', '--title',
              help     = 'Plot title.')
@click.option( '--normalize',
               is_flag = True,
               help    = 'Normalize by total metagene site coverage')
@click.option( '--dump', "-d",
               type    = click.Path(),
               help    = 'Dump the data to csv file')
def metagene(ribo, 
             site, 
             out, 
             experiments, 
             title,
             normalize,
             lowerlength, 
             upperlength, 
             dump ):
    """Generates metagene plots.
    
    The x-axis is the relative nucleotide positions
    and the start /stop site is at the center (origin at 0).
    The y-axis is the raw or normalized coverage.

    At most 7 experiments can be provided for a single plot.

    For a given start or stop site, the coverage around
    the site of interest is plotted. 

    The supported output extensions are pdf and png.

    If a length range is not provided, the minimum and the maximum
    ranges from the ribo file are read and used as the range.
    The values in the length range are aggregated to generate
    the metagene plot.

    \b
     
    Examples:  
     1) Get coverage data, of WT, for read length 21 in zipped bedgraph format

       .. code:: bash 

               ribopy plot metagene -s start -o hela_1.pdf project.ribo HeLa_1


     2) Two experiments in one plot, stop site, with normalization

       .. code:: bash 

               ribopy plot metagene -s start --normalize -o hela_1_2.pdf project.ribo HeLa_1 HeLa_2
               
     3) Plot start site usiong footprints of length 20,21,22
     Dump the data to out.csv file.

       .. code:: bash 

               ribopy plot metagene -s start \\
                 --lowerlength 20 upperlength 22 \\
                 -o hela_1.pdf \\ 
                 -d out.csv \\
                 project.ribo HeLa_1             
    """

    return plot_metagene_wrapper( 
                           ribo_file        = ribo, 
                           output_file      = out,
                           site_type        = site, 
                           range_lower      = lowerlength, 
                           range_upper      = upperlength,
                           experiment_list  = experiments,
                           title            = title,
                           normalize        = normalize,
                           dump_to_file     = dump)




@plot.command()
@click.argument("ribo",      type = click.Path(exists=True))
@click.argument("experiments", nargs = -1,)
@click.option('-r', '--region',
              type     = click.Choice(EXTENDED_REGION_names),
              required = True,
              help     = 'Region type.')
@click.option('-o', '--out',
              required = True,
              help     = 'Output file in bed format')
@click.option('-t', '--title',
              help     = 'Plot title.')
@click.option( '--normalize',
               is_flag = True,
               help    ='Normalize by total metagene site coverage')
@click.option( '--dump', "-d",
               type    = click.Path(),
               help    = 'Dump the data to csv file')
def lengthdist(ribo, out, experiments, title,
                   normalize, region, dump ):
    """Plots the distribution of the ribosome footprint lengths.

    The x-axis is the length of the protected ribosome footprints.
    The y-axis is the raw or normalized frequecies.

    At most 7 experiments can be provided for a single plot.

    Pdf and png output formats are supported.
    If "dump" option is provided, the data is written 
    to the provided file path.

    If the frequencies are normalized using the "--normalize" option, 
    the y-axis becomes the percentages of the frequencies. 

    \b
     
    Examples:  
     1) Plot CDS length distribution of exp_1 and exp_2 
     and normalize the frequencies.

       .. code:: bash 

               ribopy plot lengthdist  \\
                  -o multiple_dist.pdf \\
                  -r CDS --normalize \\
                  project.ribo exp_1 exp_2

     2)Plot only main_exp and write the data to out.csv.

       .. code:: bash 

               ribopy plot lengthdist \\
                   -d out.csv \\
                   -o main_exp.pdf \\
                   -r CDS \\
                   project.ribo  main_exp

    """

    return plot_lengthdist_wrapper( 
                           ribo_file        = ribo, 
                           output_file      = out,
                           region_type      = region, 
                           experiment_list  = experiments,
                           title            = title,
                           normalize        = normalize,
                           dump_to_file     = dump)


@plot.command()
@click.argument("ribo",      type = click.Path(exists=True))
@click.argument("experiments", nargs = -1,)
@click.option('-o', '--out',
              required = True,
              help     = 'Output file')
@click.option('--lowerlength',
              required = False ,
              type     = click.INT,
              help     = 'Lower read length')
@click.option('--upperlength',
              required = False ,
              type     = click.INT,
              help     ='Upper read length')
@click.option('-t', '--title',
              help     = 'Plot title.')
@click.option( '--horizontal',
               is_flag = True,
               help    ='Draw bars horizontally.')              
@click.option( '--dump', "-d",
               type    = click.Path(),
               help    = 'Dump the data to csv file')
def regioncounts(ribo, out, experiments, 
                 upperlength, lowerlength, title, horizontal, dump):
    """Generates barplots of the percentages of the UTR5, CDS and UTR3 counts.

    The raw counts are saved to a csv file if "--dump" option is provided.

    If a length range is not provided, then the range is determined
    using the mininum and the maximum read lengths in the ribo file.

    \b
     
    Examples:  
     1) Plot the region counts for the experiment names sample
     for lengths from 29 to 31

      .. code:: bash 

         ribopy plot regioncounts --lowerlength 29 --upperlength 31 \\
             --out sample.region_counts.pdf test.ribo sample

     2) Plot region counts for the two experiments WT and DrugA
     for all lengths combined

      .. code:: bash 

         ribopy plot regioncounts --out wt_anddrug.pdf other.ribo WT DrugA 
    """

    return plot_region_counts_wrapper(
                        ribo_file        = ribo, 
                        experiment_list  = experiments, 
                        range_lower      = lowerlength, 
                        range_upper      = upperlength, 
                        title            = title,
                        output_file      = out, 
                        dump_to_file     = dump,
                        horizontal       = horizontal)
