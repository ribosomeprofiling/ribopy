from .main import *
from ..rnaseq import *

@cli.group()
def rnaseq():
    """
    Display, set or delete RNA-Seq data 
    """
    pass


@rnaseq.command()
@click.argument('ribo', type = click.Path( ))
@click.option('-n', '--name',  
               help     = "experiment name", 
               type     = click.STRING , 
               required = True)
@click.option('-a', '--alignment',  
               help     = "RNASeq alignments in bed or bam format.", 
               type     = click.Path(exists = True) , 
               required = False)
@click.option('-c', '--counts',  
               help     = "Transcript Expression File", 
               type     = click.Path(exists = True) , 
               required = False)
@click.option('-f', '--format',  
               help     = "RNA-Seq alignment format", 
               type     = click.Choice(["bed", "bam"]) ,
               default  = "bed",
               required = False)
@click.option('--sep',  
               help         = "Column Separator for counts file", 
               type         = click.STRING ,
               default      = "\t",
               show_default = True,  
               required     = False)
@click.option( '--force',
               is_flag = True,
               help    = 'Set RNA-Seq without prompting user.')
def set(ribo, name, alignment, counts, sep, format, force):
    """
    Store the transcript expression data of an experiment

    ?? MISSING DOCUMENTATION ??

    """

    set_rnaseq_wrapper(ribo_file     = ribo, 
                       name          = name, 
                       rnaseq_file   = alignment,
                       rnaseq_counts = counts,
                       sep           = sep,
                       format        = format,
                       force         = force)


@rnaseq.command()
@click.argument('ribo', type = click.Path( ))
@click.option('--name',  
               help     = "experiment name", 
               type     = click.STRING , 
               required = False)
@click.option('--out',  
               help     = "Output File", 
               type     = click.Path() , 
               required = False)
@click.option('--sep',  
               help     = "Column Separator: Default is tab", 
               type     = click.STRING ,
               default  = "\t",  
               required = False)
def get(ribo, name, out, sep):
    """
    Get transcript expression data of a given experiment

    If no output parameter is provided, the results are printed to standard output.

    Transcript expression is reported in two columns 
    where the first column corresponds to transcript names and
    the second column corresponds to transcript expression.

    \b
    Examples
    --------

    Save the transcript expression in a tab separated file in gzipped form.

    1) ribopy rnaseq get --name WT --out transcript_exp.tsv.gz test.ribo

    Print the transcript expression on the screen

    2) ribopy rnaseq get --name WT test.ribo
    """

    get_rnaseq_wrapper(ribo_file = ribo, 
                       name      = name,
                       output    = out,
                       sep       = sep)



@rnaseq.command()
@click.argument('ribo', type = click.Path( ))
@click.option('--name',  
               help     = "experiment name", 
               type     = click.STRING , 
               required = True)
@click.option( '--force',
               is_flag = True,
               help    = 'Delete RNA-Seq without prompting user.')
def delete(ribo, name, force):
    """
    Delete RNA-Seq data of a particular experiment

    \b
    Example
    -------

    ribopy rnaseq delete --name WT test.ribo 
    """

    delete_rnaseq_wrapper(ribo_file = ribo, 
                          name      = name,
                          force     = force)
