# -*- coding: utf-8 -*-
from collections import OrderedDict
import os
from io import IOBase, StringIO

import h5py
import yaml
import numpy as np
import pandas as pd

from .settings import *
from .io.file import flex_open, flex_out, get_alignment_file_handle
from .core.get_gadgets import *
from .core.exceptions import *
from .core.verify import prompt_user, make_cli_function,\
                         check_experiment_list_in_ribo_handle
from .core.coverage import find_coverage
from .core.region_counts import find_region_counts
from ._version import (__format_version__,
                      __version__) 

###########################################################

col_types = {0: str,
             UTR5_name:RNASEQ_DT,
             UTR5_JUNCTION_name:RNASEQ_DT, 
             CDS_name:RNASEQ_DT, 
             UTR3_JUNCTION_name:RNASEQ_DT, 
             UTR3_name:RNASEQ_DT}

table_col_number_err_msg = "Invalid table!\n" + \
                           "Transcript abundance table " + \
                           "should have TRANSCRIPT_name and 5 columns for:\n" +\
                           "  UTR5 UTR5_junction " + \
                           "CDS UTR3_junction UTR3"

###########################################################
def _preprocess_rnaseq_table(region_counts, name, sep = "\t"):
    """
    Performs checks and coversion on the rnaseq table if necessary.
     
    Checks if the rnaseq table has sufficiently many columns.
    Raises error if not.
    
    If the index is two dimensional (experiment_name, transcript_name)
    then it makes it one dimensional (experiment_name).
    
    It returns a ps.DataFrame 
    """
    if type(region_counts) == pd.DataFrame:
        first_index = region_counts.index[0]
        if len(first_index) == 2:
            region_counts = region_counts.loc[name]
            
        if len(region_counts.columns) != 5:
            raise RiboBaseError(table_col_number_err_msg)
        return region_counts
    
    # Region counts must be coming from a file path or stream
    # from this point on.
    
    with flex_open(region_counts) as region_stream:
        rnaseq_lines = region_stream.readlines()
        
    upper_first_line = list( map(lambda x: x.upper(), rnaseq_lines[0].strip().split(sep)) )
    
    row_offset = 0
    col_offset = len(rnaseq_lines[1].split(sep)) - 6
    
    if CDS_name.upper() in upper_first_line and \
       UTR5_name.upper() in upper_first_line and \
       UTR3_name.upper() in upper_first_line:
       row_offset = 1
       
    if col_offset < 0:
        raise RiboBaseError(table_col_number_err_msg)
                            
    selected_entries = []
    
    for row in rnaseq_lines[row_offset:]:
        contents = row.strip().split(sep)
        selected_entries.append( sep.join(contents[col_offset:]) )
        
    fixed_stream = StringIO("\n".join(selected_entries))
        
    return pd.read_csv( fixed_stream,
                        sep       = sep,
                        names     = EXTENDED_REGION_names,
                        header    = None,
                        index_col = 0,
                        dtype     = col_types) 
    
def _get_rnaseq_from_df( region_counts, reference_list, name, sep = "\t" ):
    
    rnaseq_df           = _preprocess_rnaseq_table(
                                       region_counts = region_counts, 
                                       name          = name, 
                                       sep           = sep)
    user_reference_list = list(rnaseq_df.index)
                                         
    np_initial    = np.zeros(shape = (len(reference_list), 5), 
                             dtype = RNASEQ_DT)                     
    region_counts = pd.DataFrame(np_initial, 
                                     index   = reference_list,
                                     columns = EXTENDED_REGION_names, 
                                     dtype   = RNASEQ_DT)
                                     
                                    
    for reference in user_reference_list:

        if reference not in reference_list:
            raise ReferecenDoesntExist("{} doesn't exist ".format(reference)+\
                     "as reference in this ribo.")
        region_counts.loc[reference] = rnaseq_df.loc[reference]
  
    return region_counts 
        

@check_experiment_list_in_ribo_handle
def set_rnaseq(ribo_handle, name, rnaseq_reads, format,
               rnaseq_counts = None, sep = "\t"):
    """
    Set the rna-seq data counts.
    
    RNA-Seq data can be provided in one of the two ways:
       1) Alignments
       2) Counts by the regions 
    
       Alignments:
          Alignments can be provided in bam or bed format.
          If RNA-Seq data is paired end, users should provide alignments coming
          from the first roun of the sequencing (read 1, R1).
          
       Counts by the regions:
          If transcript abundance is already determined, users can provide the
          counts table in a tab (or comma) separated file format.
          For each region (UTR5, UTR5_junction, CDS, UTR3_junction and UTR3),
          abundance values must be provided in a separate column.
          Thus this table must be of the following form:
          
             ===============  ====  =============  ===  =============  ====
             TRANSCRIPT_NAME  UTR5  UTR5_junction  CDS  UTR3_junction  UTR3
             ---------------  ----  -------------  ---  -------------  ----
             GAPDH             2.5           7.2    78             4    1.2
             ===============  ====  =============  ===  =============  ====
          
    ?? This documentation is incomplete ??
            
    Parameters
    ----------
    ribo_handle: h5py.File
        Open ribo file handle

    name: str
        Name of the experiment

    rnaseq_df: pd.DataFrame
        Index is the transcript names
        The first (and only) column has transcript expression values.

    Returns
    -------
    actual_rnaseq_df[0]: np.array
       Transcript expression array
    """
    
    reference_names   = get_reference_names(ribo_handle)
    reference_lengths = get_reference_lengths(ribo_handle)
    annotation        = get_region_boundaries(ribo_handle)
    
    left_span, right_span = get_left_right_spans(ribo_handle)

    if rnaseq_reads:
        with get_alignment_file_handle(alignment_file = rnaseq_reads, 
                                  alignment_format = format ) as rnaseq_stream:
            
            rnaseq_coverage = find_coverage(rnaseq_stream, 
                                      ref_names   = reference_names, 
                                      ref_lengths = reference_lengths)
                                      
            region_counts = find_region_counts( coverage   = rnaseq_coverage, 
                                                annotation = annotation, 
                                                left_span  = left_span, 
                                                right_span = right_span )

    elif rnaseq_counts:
        region_counts = _get_rnaseq_from_df(rnaseq_counts, 
                                            reference_names,name, sep = sep)
    else:
        raise( RiboBaseError("set_rnaseq function needs either "
                             "rnaseq_ counts as dataframe "
                             "or rnaseq reads as input") )
                
    experiment_handle = ribo_handle[EXPERIMENTS_name][name]
    rnaseq_datagroup  = experiment_handle.require_group(RNASEQ_name)
    rnaseq_dataset    = rnaseq_datagroup.get(RNASEQ_name, None)
    
    if rnaseq_dataset:
        rnaseq_dataset[...] = region_counts
    else:
        rnaseq_datagroup.create_dataset(
                                RNASEQ_name, 
                                shape       = (len(reference_names), 5),
                                data        = region_counts,
                                dtype       = RNASEQ_DT,
                                compression = DEFAULT_COMPRESSION, 
                                fletcher32  = DEFAULT_FLETCHER32)
                                


@make_cli_function
def set_rnaseq_wrapper(ribo_file, name, rnaseq_file, rnaseq_counts = None,
                       format = "bed", sep = "\t", force = True):
    """
    Wrapper for set_metadata.
    ?? Doc Needs Revision ??

    It makes sure that the provided transcript names
    are actually in the reference transcripts.

    See set_rnaseq for details,
    """
    
    if not (rnaseq_file or rnaseq_counts):
        raise( RiboBaseError("Setting RNA-Seq data requires either "
                             "rnaseq_counts "
                             "or rnaseq reads as input") )
 
    with h5py.File(ribo_file, "r+") as ribo_handle:

        experiment_handle = ribo_handle[EXPERIMENTS_name][name]
        prompt_message    = "This will overwrite existing RNA-Seq data." 
        attr_exists       = \
           experiment_handle.get(RNASEQ_name , None) \
              and \
           experiment_handle[RNASEQ_name].get(RNASEQ_name, name)

        prompt_user( message     = prompt_message, 
                     attr_exists = attr_exists, 
                     force       = force)
                    
        set_rnaseq(ribo_handle   = ribo_handle, 
                   name          = name, 
                   rnaseq_reads  = rnaseq_file, 
                   format        = format,
                   rnaseq_counts = rnaseq_counts, 
                   sep           = sep)


def _get_single_rna_seq(ribo_handle, name):
    """
    Helper function for get_rnaseq

    It return rnaseq data of a single experiment if exists
    """

    if not rnaseq_exists(ribo_handle, name):
        raise NORNASEQ("{} doesn't have RNA_seq data".format(name))

    rnaseq_np = ribo_handle[EXPERIMENTS_name]\
                    [name][RNASEQ_name][RNASEQ_name][...]
    ref_names = get_reference_names(ribo_handle)
    
    rna_index_pre = list(zip(  [name] * len(ref_names), ref_names  )) 
    rna_index     = pd.MultiIndex.from_tuples( rna_index_pre  ) 
    rnaseq_df = pd.DataFrame( rnaseq_np , 
                              index   = rna_index , 
                              columns = EXTENDED_REGION_names)

    return rnaseq_df


def _get_all_rna_seq(ribo_handle):
    """
    Helper function for get_rnaseq

    It return rnaseq data of a all experiments in the ribo file
    """

    has_rnaseq = False

    experiment_list = get_experiment_names(ribo_handle)
    rnaseq_df_list  = list()

    for experiment in experiment_list:
        if rnaseq_exists(ribo_handle, experiment):
            has_rnaseq = True
            rnaseq_df_list.append( _get_single_rna_seq(ribo_handle, experiment) )

    if not has_rnaseq:
        return None
    else:
        all_df = pd.concat(rnaseq_df_list, axis = 0)
        return all_df

@check_experiment_list_in_ribo_handle
def get_rnaseq(ribo_handle, name = None):
    """
    Returns rnaseq data if it exists

    Parameters
    ----------
    ribo_handle: h5py.File
        Open ribo file handle

    name: str
        Name of the experiment

    Returns
    -------
    rnaseq_df: pd.DataFrame 
        transcript expression table.
        Indices are transcript names
    """

    if not name:
        rnaseq_df =  _get_all_rna_seq(ribo_handle)
    else:
        rnaseq_df = _get_single_rna_seq(ribo_handle, name)
        
    rnaseq_df.index.set_names([EXPERIMENTS_name, REFERENCE_name], 
                               inplace = True)
    
    return rnaseq_df



@make_cli_function
def get_rnaseq_wrapper(ribo_file, name, output, sep = "\t"):
    """
    Wrapper for get_rnaseq
    """

    with h5py.File(ribo_file , "r") as ribo_handle:

        try:
            rnaseq_df = get_rnaseq(ribo_handle, name)
        except NORNASEQ as e:
            print(e)
            exit(1)

        out_str = rnaseq_df.to_csv(sep = sep, header = True)
     
        with flex_out(output) as output_stream:
            print(out_str, file = output_stream)



@check_experiment_list_in_ribo_handle
def delete_rnaseq(ribo_handle, name):
    """
    Returns rnaseq profile of the data if it exists

    Parameters
    ----------
    ribo_handle: h5py.File
        Open ribo file handle

    name: str
        Name of the experiment

    Returns
    -------
    rnaseq_df: pd.DataFrame 
        transcript expression table.
        Indices are transcript names
    """

    if not rnaseq_exists(ribo_handle, name):
        raise NORNASEQ("{} doesn't have RNA_seq data".format(name))

    del ribo_handle[EXPERIMENTS_name][name][RNASEQ_name]


@make_cli_function
def delete_rnaseq_wrapper(ribo_file, name, force = True):
    """
    Wrapper for delete_rnaseq
    """

    prompt_message    = "This will delete existing RNA-Seq data."

    with h5py.File(ribo_file , "r+") as ribo_handle:
        attr_exists = rnaseq_exists(ribo_handle, name)
        prompt_user( message     = prompt_message, 
                     attr_exists = attr_exists, 
                     force       = force)
        delete_rnaseq(ribo_handle, name)
