from .main import *
from ..metadata import *

@cli.group()
def metadata():
    """
    Display, set or delete user-defined metadata

    If no name is given, the metadata of the ribo file is 
    set, displayed or deleted.
    """
    pass


@metadata.command()
@click.argument('ribo', type = click.Path( ))
@click.option('--name',  
               help     = "experiment name", 
               type     = click.STRING , 
               required = False)
@click.option('--meta',  
               help     = "Metadata File", 
               type     = click.Path(exists = True) , 
               required = True)
@click.option( '--force',
               is_flag = True,
               help    = 'Set metadata without prompting user.')
def set(ribo, name, meta, force):
    """Stores the metadata in the meta file.

    If no name is given, the metadata of the ribo file is set.

    If name is provided, the metadata of the corresponding experiment
    is set.

    The metadata must be in yaml format.

    \b
    Example Ribo Metadata File Contents:
          
        =================  ===================
        pipeline_name:     RiboFlow
        pipeline_version:  v1.0.2
        project:           elongation blocker
        =================  ===================


    \b
    Example Library Metadata File Contents:
        
       ===========   ========== 
        Cell_Line:   Human ESC
        Treatment:   Drug A   
        Enzyme:      RNASEI
       ===========   ==========      

    """

    set_metadata_wrapper(ribo_file  = ribo , 
                         name       = name , 
                         meta_file  = meta,
                         force      = force)


@metadata.command()
@click.argument('ribo', type = click.Path( ))
@click.option('--name',  
               help     = "experiment name", 
               type     = click.STRING , 
               required = False)
def get(ribo, name):
    """Displays user-defined metadata of the ribo file or experiment

    If no name is provided, metadata of the ribo file is displayed.

    If name is given, metadata of the corresponding experiment is displayed.
    """

    get_metadata_wrapper(ribo_file = ribo , 
                         name      = name)


@metadata.command()
@click.argument('ribo', type = click.Path( ))
@click.option('--name',  
               help     = "experiment name", 
               type     = click.STRING , 
               required = False)
@click.option( '--force',
               is_flag = True,
               help    = 'Set metadata without prompting user.')
def delete(ribo, name, force):
    """
    Deletes user-defined metadata of the ribo file or experiment

    If no name is provided, metadata of the ribo file is deleted.

    If name is given, metadata of the corresponding experiment is deleted.
    """

    delete_metadata_wrapper(ribo_file = ribo , 
                            name      = name,
                            force     = force)
