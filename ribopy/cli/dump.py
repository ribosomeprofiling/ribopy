# -*- coding: utf-8 -*-

import click

from .main import *

from ..dump import *

##################################

@cli.group()
def dump():
    """Dump selected parts of ribo files to particular formats"""
    pass


@dump.command()
@click.argument("ribo", type = click.Path(exists = True))
@click.option('-s', '--site',
              type      = click.Choice(['start', 'stop']),
              required  = True,
              help      = 'Site type.')
@click.option('-o', '--out',
              help      = 'Output file in bed format')
@click.option('-e', '--experiment',
              help      = 'Name of the experiment')
@click.option( '-l', '--lowerlength',
               type     = click.INT,
               help     = 'Minimum read length to be taken')
@click.option( '-u', '--upperlength',
               type     = click.INT,
               help     = 'Maximum read length to be taken')
@click.option( '--sumlengths',
                is_flag = True,
                help    = 'Sum accross lengths')
@click.option( '--nosumtrans',
               is_flag  = True,
               help     = "Do NOT aggregate values"
                          " accross transcripts")
def metagene(ribo, 
             site, 
             out, 
             experiment, 
             lowerlength, upperlength, 
             sumlengths, 
             nosumtrans ):
    """Dumps metagene data 
    to a csv file or standard output.

    Metagene data is obtained from
    coverage around start / stop site in a pre-defined radius.


    Examples:
     Below are some examples for different scenarios.
     
     1) Get metagene data around START site for the experiment WT.
        Report results for read lengths from 28 to 32. 
        Agrregate data by summing accros read lengths.
        Also, values are summed across transcripts.
        Save the results in start.csv.
          
            \b 
            
            .. code:: bash      
             
             ribopy dump metagene \\
                --experiment WT \\
                --site start \\
                --out start.csv \\
                --lowerlength 28 --upperlength 32 \\
                --sumlengths \\
                sample.ribo
                
     2) Get metagene data around STOP site for the experiment Treatment_1.
        Report results for read length 31. 
        Values are summed across transcripts.
        Print the results on the standard output.
          
            \b 
            
            .. code:: bash      
             
             ribopy dump metagene \\
                --experiment WT \\
                --site stop \\
                --lowerlength 31 --upperlength 31 \\
                sample.ribo
                
     3) Get metagene data around START site for the experiment WT.
        Report results for read lengths from 30 to 32. 
        Report results for each read length.
        Also, values are summed across transcripts.
        Save the results in start.csv.
          
            \b 
            
            .. code:: bash      
             
             ribopy dump metagene \\
                --experiment WT \\
                --site start \\
                --out start.csv \\
                --lowerlength 30 --upperlength 32 \\
                sample.ribo
                
    """

    experiment_list = []
    if experiment:
        experiment_list = [experiment]

    return dump_metagene_wrapper( 
                           ribo_file       = ribo, 
                           output_file     = out,
                           site_type       = site, 
                           sum_lengths     = sumlengths, 
                           sum_references  = (not nosumtrans),
                           range_lower     = lowerlength, 
                           range_upper     = upperlength,
                           experiment_list = experiment_list)




###########################################################

@dump.command()
@click.argument("ribo",   type = click.Path(exists=True))
@click.option('-r', '--region',
              type      = click.Choice(['UTR5', 'UTR5_junction',
                                        'CDS',
                                        'UTR3', 'UTR3_junction']),
              required  = True,
              help      = "Site type.")
@click.option('-o', '--out',
              help      = 'Output file in csv format')
@click.option('-e', '--experiment',
              help      = 'Name of the experiment')
@click.option( '--lowerlength',
               type     = click.INT,
               help     = 'Minimum read length to be taken')
@click.option( '--upperlength',
               type     = click.INT,
               help     = 'Maximum read length to be taken')
@click.option( '--sumlengths',
                is_flag = True,
                help    = 'Sum accross lengths')
@click.option( '--sumtrans',
               is_flag  = True,
               help     = 'Sum accross transcripts')
def region_counts(ribo, 
                  region, 
                  out, 
                  experiment,
                  lowerlength, 
                  upperlength, 
                  sumlengths, 
                  sumtrans ):
    """
    Dumps a given region in csv format to a file or standard output.
    
    Examples:
     Below are some examples for different scenarios.
     
     1) Get number of reads mapping to the coding sequence (CDS) 
        for the experiment WT.
        Report results for read lengths from 28 to 32. 
        Agrregate data by summing accros read lengths.
        Also, values are summed across transcripts.
        Save the results in cds.csv.
          
            \b 
            
            .. code:: bash      
             
             ribopy dump region-counts \\
                --experiment WT \\
                --region CDS \\
                --out cds.csv \\
                --lowerlength 28 --upperlength 32 \\
                --sumlengths \\
                --sumtrans \\
                sample.ribo
                
     2) Get number of reads mapping to UTR3 
        for the experiment Treatment.
        Report results for read lengths from 30 to 32. 
        Agrregate data by summing accros read lengths.
        CDS occupancy is reported for each transcript.
        Save the results in cds.csv.
          
            \b 
            
            .. code:: bash      
             
             ribopy dump region-counts \\
                --experiment Treatment \\
                --region UTR3 \\
                --out cds.csv \\
                --lowerlength 30 --upperlength 32 \\
                --sumlengths \\
                sample.ribo
    """

    experiment_list = []
    if experiment:
        experiment_list = [experiment]

    return dump_region_wrapper(
                        ribo_file       = ribo, 
                        output_file     = out, 
                        region_name     = region,
                        sum_lengths     = sumlengths, 
                        sum_references  = sumtrans, 
                        range_lower     = lowerlength, 
                        range_upper     = upperlength,
                        experiment_list = experiment_list)


######################################################

@dump.command()
@click.argument("ribo", type = click.Path(exists = True))
@click.option('-o', '--out',
              help = 'Output file in bed format')
def annotation(ribo, out):
    """Gives UTR5, CDS and UTR3 annotation in bed format.

    Annotation is written to standard output if no
    output file is given.

    Examples:
     1) Store the annotation in a bed file 
     

     .. code:: bash
         
         ribopy dump annotation -o regions.bed sample.ribo
         
     2) Print annotation to standard output 
     
     .. code:: bash
         
        ribopy dump annotation sample.ribo
         
     
    """
    return dump_annotation_wrapper(ribo, out)


@dump.command()
@click.argument("ribo", type = click.Path(exists=True))
@click.option('-o', '--out',
              help     = 'Output file',
              type     = click.Path(False),
              required = False)
@click.option('-s', '--sep',
              help            = 'Column separator, default is \\t (tab)',
              type            = click.STRING,
              default         = "\t",
              required        = False)
def reference_lengths(ribo, out, sep):
    """Gives transcript names and their lengths.

    Annotation is written to standard output if no
    output file is given.

    The first column corresponds to transcript names 
    and the second column corresponds to transcript
    lengths.

    Examples:
     Print the output to the terminal

     .. code:: bash
     
         ribopy dump reference-lengths sample.ribo

     Save the output, in gzipped form, in lengths.csv.gz
     and use , to separate columns

     .. code:: bash

         ribopy dump reference-lengths -o lengths.csv.gz --sep "," sample.ribo
    """
    return dump_reference_lengths_wrapper(ribo, out, sep = sep)

#########################################################

@dump.command()
@click.argument("ribo",    type = click.Path(exists=True))
@click.argument("experiment", type = click.STRING)
@click.option( '-o', '--out',
               help     = 'Output file in csv format')
@click.option( '--lowerlength',
               type     = click.INT,
               help     = 'Minimum read length to take')
@click.option( '--upperlength',
               type     = click.INT,
               help     = 'Maximum read length to take')
@click.option( '--format',
               type     = click.Choice(['bg', 'tsv']),
               help     = 'Output file format')
def coverage(ribo, 
             out, 
             experiment,
             lowerlength, upperlength, 
             format):
   """
   Prints the coverage of each transcript.

   This command prints the coverage of each transcript, 
   at nucleotide resolution, for a given read length range.
   For a single read length, set 
   'lowerlength' equal to 'upperlength'.
   The values for the given range are aggregated by summation.

   Append ".gz" to the output file name to have the 
   output in zipped form. 

   File format:     
    options: bg, tsv

    bg: Bedgraph
     The columns of the bedgraph file are of the form
     transcript_name location_start location_end coverage_value
     Only the nonzero values are reported in the bedgraph file.
     Bedgraph is 0-based and location_start is inclusive and
     location_end is exclusive. 

     See https://genome.ucsc.edu/goldenPath/help/bedgraph.html
     for a detailed description of bedgraph file format.

    tsv: Tab separated File
      In this format, coverage values for all nucleotide positions,
      regardless of their value ( 0 or not ) are reported.

      The columns of the file are of the form
      transcript_name comma_separated_coverage 

   \b
    
   Examples:  
    1) Get coverage data, of WT, for read length 21 in zipped bedgraph format

    \b
    
    .. code:: bash 
       
               ribopy --lowerlength 21 --upperlength 21 \\
                   -o coverage.bg.gz \\
                   --fromat bg sample.ribo WT


    2) Get coverage data, of treatment_1, for range 26 to 30 in tsv format

    \b
    
    .. code:: bash 
       
           ribopy --lowerlength 26 --upperlength 30 
                 -o coverage.tsv  
                 --format tsv 
                 sample.ribo treatment_1
   """

   dump_coverage_wrapper(ribo_file       = ribo, 
                         output_file     = out, 
                         experiment_name = experiment,
                         range_lower     = lowerlength, 
                         range_upper     = upperlength,
                         file_format     = format)
