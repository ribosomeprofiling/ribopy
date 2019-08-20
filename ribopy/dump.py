# -*- coding: utf-8 -*-

from functools import partial, wraps
from collections import OrderedDict

import h5py
import pandas as pd
import numpy as np

from .settings import *
from .core.get_gadgets import *
from .core.aggregate_by_length import aggregate_region_counts,\
                                      aggregate_start_site,\
                                      aggregate_stop_site

from .core.verify import verify_extended_region_name,\
                         verify_site_type

from .core.exceptions import *

from .io.file import flex_out
from .io.dump_coverage import print_coverage

from .core.verify import make_cli_function, \
                         check_experiment_list_in_ribo_handle, \
                         check_length_range

######################################################################                         

COLUMN_INDECES = dict( zip(EXTENDED_REGION_names, 
                           tuple(range(0, len(EXTENDED_REGION_names)))))

#######################################################################

#######################################################################

@check_experiment_list_in_ribo_handle
@check_length_range
def dump_metagene( ribo_handle, *,
                   site_type, 
                   sum_lengths, 
                   sum_references,
                   range_lower     = 0, 
                   range_upper     = 0,
                   experiment_list = [],
                   alias           = None ):
    """
    **TODO:**
    Complete documentation
    """

    length_min, length_max = get_read_length_range(ribo_handle)

    if not range_lower:
        range_lower = length_min
    if not range_upper:
        range_upper = length_max

    verify_site_type(site_type)

    if not experiment_list:
        experiment_list = get_experiment_names(ribo_handle)

    if site_type == REF_DG_START_SITE_COV:
        get_site = aggregate_start_site
    elif site_type == REF_DG_STOP_SITE_COV:
        get_site = aggregate_stop_site
    else:
        raise ValueError("Wring site type "+ site_type +\
                         "It must be one of {} or {}".format(\
                           REF_DG_START_SITE_COV, REF_DG_STOP_SITE_COV ))

    metagene_data = get_site(ribo_handle     = ribo_handle, 
                             range_lower     = range_lower, 
                             range_upper     = range_upper,
                             sum_values      = sum_lengths,
                             experiment_list = experiment_list)

    metagene_radius = get_metagene_radius(ribo_handle)

    # Note that the second value is excluded
    # so we add +1 tp include the +metagene_radius position
    coverage_interval = list(range( -1 * metagene_radius, 
                                     metagene_radius + 1 ) )
    reference_names = list(get_reference_names(ribo_handle, alias))

    if sum_lengths:
        column_labels = coverage_interval
        create_df_pre = partial( pd.DataFrame, columns = column_labels)

        if sum_references:
            create_df = (lambda x: pd.DataFrame(
                               create_df_pre(x).sum()   ) )
        else:
            create_df = partial( pd.DataFrame, columns = column_labels,
                                 index = reference_names)
    else:
        column_labels = [DF_READLENGTH]  + coverage_interval
        create_df_pre = partial( pd.DataFrame, columns = column_labels)

        if sum_references:
            create_df = (lambda x: create_df_pre(x).\
                            groupby([DF_READLENGTH]).sum()  )
        else:
            reference_name_repetition_time = (range_upper - range_lower) + 1 
            reference_labels = \
                reference_names * reference_name_repetition_time
            create_df = partial( pd.DataFrame, columns = column_labels,
                                     index = reference_labels)


  
    metagene_dataframes = list(map(create_df , metagene_data))

    if sum_lengths and sum_references:
        metagene_dataframes = list(map(lambda x: x.T, 
                                 metagene_dataframes) )

    dataframes_dict = OrderedDict( 
                          tuple(zip(experiment_list, 
                                     metagene_dataframes)) )

    return dataframes_dict




def merge_metagene_dataframes(df_dict, 
                              sum_references, 
                              sum_lengths):

    """
    Incomplete documentation
    """

    experiment_list = tuple(df_dict.keys())

    if sum_references:
        names = [DF_EXPERIMENT_NAME]
    else:
        names = [DF_EXPERIMENT_NAME, DF_TRANSCRIPT]

    # We try to avoid the redundant index colum
    # if references and lengths are summed
    if sum_references and sum_lengths:
        merged_df = pd.concat( 
                       tuple(df_dict.values()), 
                       keys  = experiment_list,
                       names = names,
                       join  = "outer").reset_index(level=1, drop=True)

    else:
        merged_df = pd.concat( 
                       tuple(df_dict.values()), 
                       keys  = experiment_list,
                       names = names,
                       join  = "inner")
                       
    # If lengths and references are not summed,
    # length will stay in the columns not indices.
    # So we move lengths to the index part
    # by reindexing the dataframe.                     
    if (not sum_lengths) and (not sum_references):
        merged_df =  merged_df.set_index([merged_df.index, DF_READLENGTH])
        
    return merged_df

@make_cli_function
def dump_metagene_wrapper( ribo_file, 
                           output_file, 
                           site_type, 
                           sum_lengths, 
                           sum_references,
                           range_lower     = 0, 
                           range_upper     = 0,
                           experiment_list = [] ):
    """
    TODO:
    Complete documentation
    """

    if site_type.lower() == "start":
        site_type = REF_DG_START_SITE_COV
    elif site_type.lower() == "stop":
        site_type = REF_DG_STOP_SITE_COV
    else:
        raise ValueError("Wring site type:", site_type)


    ribo_handle = h5py.File(ribo_file, "r")

    try:
        df_dict = dump_metagene( ribo_handle     = ribo_handle,
                                 site_type       = site_type, 
                                 sum_lengths     = sum_lengths, 
                                 sum_references  = sum_references,
                                 range_lower     = range_lower, 
                                 range_upper     = range_upper,
                                 experiment_list = experiment_list ) 

        site_coverage = merge_metagene_dataframes(df_dict, 
                                 sum_references, sum_lengths)
    except ValueError as e:
        print(e)
        exit(1)
    finally:
        ribo_handle.close()

    if output_file:
        site_coverage.to_csv(output_file, sep = DF_CSV_SEPARATOR)
    else:
        print(site_coverage.to_csv( sep = DF_CSV_SEPARATOR) )


#########################################################################

@check_experiment_list_in_ribo_handle
@check_length_range
def dump_region(ribo_handle, *,
                region_name,
                sum_lengths, 
                sum_references, 
                range_lower      = 0 , 
                range_upper      = 0 ,
                experiment_list  = [],
                alias            = None):

    """
    Dumps a selection region counts to a dataframe.
    Available regions are:
       UTR5, UTR5_junction, CDS, UTR3_junction, UTR3

    For each experiment, there will be a column with the experiment name.

    Originally, we keep counts for all transcripts and all lengths
    The counts can be summed across transcripts or across lengths
    by setting the flags sum_lengths or sum)references.

    If the counts are NOT summed across lengths, then
    the first column will be for the corresponding read lengths

    If the sums are not across transcripts (references), then
    the index column will indicate the corresponding transcript (reference)
    name.

    Parameters
    ----------

    ribo_handle : h5py.File
        A pointer to open hdf5 file

    region_name : str
            Type of the region. 
            Choices are: UTR5, UTR5_junction, CDS, UTR3_junction, UTR3

    sum_lengths : Boolean
            sum the counts across the lengths

    sum_references : Boolean
            sum the counts across transcripts (references)

    range_lower : int
            starting read length to be taken

    range_upper : int
            ending read length to be taken

    experiment_list : list(str)
                   the experiments to be reported
                   If empty, all experiments will be taken.

    Returns
    -------

    result_df : A pandas dataframe

    """

    length_min, length_max = get_read_length_range(ribo_handle)

    if not range_lower:
        range_lower = length_min
    if not range_upper:
        range_upper = length_max

    verify_extended_region_name(region_name)

    if experiment_list:
        actual_experiment_list = experiment_list
    else:
        actual_experiment_list = get_experiment_names(ribo_handle)
    
    region_counts = aggregate_region_counts( 
                             ribo_handle     = ribo_handle, 
                             range_lower     = range_lower, 
                             range_upper     = range_upper, 
                             sum_values      = sum_lengths,
                             experiment_list = actual_experiment_list )



    
    # Note that the first column is for the read length
    # so we need to add 1 to the column index
    offset           =  0 if sum_lengths else 1 
    collected_counts = [ c[:, COLUMN_INDECES[region_name] + offset]
                            for c in region_counts ]

    reference_names  = get_reference_names(ribo_handle, alias = alias)
    df_index         = []
    column_labels    = []

    #### If we are NOT summing the counts by read lengths
    #### Then the counts are reported per read length
    #### So we need to keep them in the first column
    if not sum_lengths:
        read_lengths  = region_counts[0][:,0]
        df_contents   = [read_lengths] + collected_counts
        column_labels = [DF_READLENGTH] + list(actual_experiment_list)
    else:
        column_labels = actual_experiment_list
        df_contents   = collected_counts


    data_table = np.array( df_contents )
    counts_df  = pd.DataFrame(np.transpose(data_table),
                             columns = column_labels)

    if sum_references:
        if sum_lengths:
            result_df = pd.DataFrame(counts_df.sum()).T
        else:
            result_df = counts_df.groupby([DF_READLENGTH]).sum()
    else:
        if sum_lengths:
            reference_name_repetition_time = 1
            reference_labels               = list(reference_names) *\
                                               reference_name_repetition_time
            column_labels                  = actual_experiment_list
            
            result_df  = pd.DataFrame(np.transpose(data_table), 
                                      index   = reference_labels,
                                      columns = column_labels).\
                                         rename_axis( DF_TRANSCRIPT)
        else:
            reference_name_repetition_time = (range_upper - range_lower) + 1
            reference_labels               = list(reference_names) *\
                                                reference_name_repetition_time
            column_labels                  = \
                 [DF_READLENGTH] + list(actual_experiment_list)
                 
            result_df  = pd.DataFrame(np.transpose(data_table), 
                                      index   = reference_labels,
                                      columns = column_labels).\
                                         rename_axis( DF_TRANSCRIPT)
            result_df = result_df.set_index([ result_df.index, DF_READLENGTH ])

    return result_df




@make_cli_function   
def dump_region_wrapper(ribo_file, 
                        output_file, 
                        region_name,
                        sum_lengths, 
                        sum_references, 
                        range_lower     = 0, 
                        range_upper     = 0,
                        experiment_list = []):

    """
    A wrapper function for the dump_region function.

    If output file is given, the resulting data frame will be saved
    in it. Else, the dataframe will be printed on the standard output
    in csv format.

    See dump_region function for details.
    """

    ribo_handle = h5py.File(ribo_file, "r")

    try:
        region_counts = dump_region( ribo_handle     = ribo_handle, 
                                     region_name     = region_name,
                                     sum_lengths     = sum_lengths, 
                                     sum_references  = sum_references, 
                                     range_lower     = range_lower , 
                                     range_upper     = range_upper ,
                                     experiment_list = experiment_list)
    except ValueError as e:
        print(e)
        exit(1)
    finally:
        ribo_handle.close()

    if output_file:
        region_counts.to_csv(output_file, sep = DF_CSV_SEPARATOR)
    else:
        print(region_counts.to_csv( sep = DF_CSV_SEPARATOR) )



##########################################################

def dump_annotation(ribo_handle):
    """
    Returns annotation of a ribo file in bed format
    in string form.

    Parameters
    ----------

    ribo_handle : h5py.File
         hdf5 handle for the ribo file


    Returns
    -------

    A string that can be output directly as a bed file.
    """

    boundaries = get_region_boundaries(ribo_handle)
    names      = get_reference_names(ribo_handle)

    bed_rows   = list()

    for ref_name, ref_boundaries in zip(names, boundaries):
        for region_name, region_boundaries in zip(REGION_names, ref_boundaries ):
            if region_boundaries[1] <= region_boundaries[0]:
                continue
            bed_entries = tuple( map( str, 
                                 [ref_name, 
                                  region_boundaries[0], region_boundaries[1],
                                  region_name, 0, "+"] ) )
            bed_rows.append( "\t".join(bed_entries) )
    return "\n".join(bed_rows)



@make_cli_function
def dump_annotation_wrapper(ribo_file, output_file = ""):
    """
    Dumps region annotation of a given ribo file.
    If the output file is not provided, 
    the result is printed on the standard output

    parmaters
    ---------
    ribo_file: ribo file path

    output_file: output file path
    """

    ribo_handle = h5py.File(ribo_file, "r")
    annotation = dump_annotation(ribo_handle)

    if output_file:
        with open(output_file, "w") as output_stream:
            print(annotation, file=output_stream)
    else:
            print(annotation)

    ribo_handle.close()

#####################################################################

def dump_reference_lengths(ribo_handle):
    """
    Dumps refrence names and their lengths in a given
    ribo file.

    The output is a pandas dataframe whose index is reference names

    parameters
    ----------
    ribo_file: h5py.File
       ribo handle

    Returns
    -------
    ref_lens_df: pandas.DataFrame
    """

    ref_lengths  = get_reference_lengths(ribo_handle)
    ref_names    = get_reference_names(ribo_handle)
    ref_lens_df  = pd.DataFrame(ref_lengths , 
                                 index   = ref_names,
                                 columns = ["length"], 
                                 dtype   = TRANSCRIPT_COVERAGE_DT) 
    return ref_lens_df


@make_cli_function
def dump_reference_lengths_wrapper(ribo_file, 
                                   output_file = "", sep = "\t"):
    """
    Dumps refrence names and their lengths in a given
    ribo file.

    If no output file is given, the result is printed to
    the standard output.

    parameters
    ----------
    ribo_file: ribo file path

    output_file: output file path
    """

    ribo_handle       = h5py.File(ribo_file, "r")
    ref_lengths_df    = dump_reference_lengths(ribo_handle)
    ref_lengths_str   = ref_lengths_df.to_csv(sep = sep, header = None)

    with flex_out(output_file) as output_stream:
        print(ref_lengths_str, file = output_stream)

    ribo_handle.close()


#####################################################################

@check_experiment_list_in_ribo_handle
@check_length_range
def dump_coverage( ribo_handle, *, 
                   experiment_name,
                   range_lower = 0, 
                   range_upper = 0,
                   alias       = None):

    """
    Returns an ordered dictionary 
    whose keys are transcript names 
    and values are coverage of the transcript.

    The coverage values in the range
    range_lower and range_upper
    are summed up.
    So if N = range_lower = range_upper,
    then the coverage values at length N
    are reported. If no range is provided,
    the reported range goes from the minimum length to maximum length
    of the ribo file.  

    Arguments
    ---------
    ribo_handle: h5py.File
       Handle to open ribo file

    experiment_name: str
       Name of the experiment 
       whose coverage is to be extracted

    range_lower: int
       Minimum read length to be taken

    range_upper: int
       Maximum read length to be taken 

    Returns
    -------
    transcript_coverage: OrderedDict( np.array )
       The keys of transcript coverage are transcript names
       The values are numpy arrays giving the coverage.
       So the value at position i corresponds to coverage on
       the ith nucleotide.
    """

    if not has_coverage_data(ribo_handle, experiment_name):
        raise ExperimentDoesntExist(
               "The experiment {} doesn't have coverage data"\
                .format(experiment_name))


    ref_lengths = get_reference_lengths(ribo_handle)
    ref_names   = get_reference_names(ribo_handle, alias)

    min_length, max_length = \
        get_read_length_range(ribo_handle)

    total_nucleotides  = sum(ref_lengths)
    total_slice_array  = np.zeros( total_nucleotides, 
                                   dtype = TRANSCRIPT_COVERAGE_DT )
    length_offset      = range_lower - min_length
    coverage_handle    = ribo_handle[EXPERIMENTS_name]\
                                    [experiment_name][REF_DG_COVERAGE]

    for this_length in range(range_lower, range_upper + 1):
        start_index = (this_length - min_length) * total_nucleotides
        end_index   = start_index + total_nucleotides  

        total_slice_array += \
            coverage_handle[REF_DG_COVERAGE][start_index : end_index]

    transcript_coverage = OrderedDict()
    start_index         = 0
    end_index           = 0

    for t_name, t_length in zip(ref_names, ref_lengths):
        end_index                   = start_index + t_length 
        transcript_coverage[t_name] = \
              total_slice_array[start_index : end_index ]
        start_index += t_length

    return transcript_coverage

########################################################################


@make_cli_function
def dump_coverage_wrapper(ribo_file, 
                          output_file, 
                          experiment_name,
                          range_lower = 0, 
                          range_upper = 0,
                          file_format = "bg",
                          separator   = ","):

    """
    This is a wrapper function for dump_coverage.

    This function arranges the output files and
    calling appropriate io functions for writing
    the output.

    See the documentation of dump_coverage for details.
    """

    ribo_handle = h5py.File(ribo_file, "r")

    if file_format not in ["bg", "tsv"]:
        print("File format can be either \"bg\" or \"tsv\".")
        exit(1)

    if not has_coverage_data(ribo_handle, experiment_name):
        exit("Error: The experiment {} does not have coverage data."\
              .format(experiment_name))

    with flex_out(output_file, "wt") as outstream:

        coverage =  dump_coverage( 
                       ribo_handle     = ribo_handle, 
                       experiment_name = experiment_name,
                       range_lower     = range_lower, 
                       range_upper     = range_upper)

        print_coverage(coverage, 
                        file_handle = outstream,
                        file_format = file_format,
                        sep         = separator)

    ribo_handle.close()

##########################################################################    
