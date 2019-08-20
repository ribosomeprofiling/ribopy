# -*- coding: utf-8 -*-
from functools import partial

from .get_gadgets import *
from ..settings import *

###########################################################

def sum_individual_dataset_list(dataset, 
                                number_of_references, 
                                lower_relative_index, 
                                upper_relative_index):
    """
    Sums the counts coming from an individual dataset
    of one experiment.

    For each reference (transcript), it grabs the counts in the
    given range and sums them and reports the sum.
    These inintial counts, being summed, are typically coming 
    from the hdf5 dataset directly.

    Both range values are inclusive.

    Parameters
    ----------

    dataset: Actual dataset of count values from hdf5

    number_of_references: number of references (transcripts)
                          in the annotation of the ribo file

    lower_relative_index: relative position of the initial
                          length relative to the minimum 
                          read length of the ribo file
                          i.e.: length_range_min - ribo_length_min

    upper_relative_index: relative position of the terminal
                          length relative to the maximum 
                          read length of the ribo file
                          i.e.: lower_relative_index + \
                                (range_max - range_min)

    Returns
    -------
    result: 2D numpy array of dimension 
            ( number_of_transcripts(references), ncols_of_dataset )

    """
    
    if len(dataset.shape) == 2:
        ncols  = dataset.shape[1]
        result = np.zeros( (number_of_references, ncols), 
                           dtype = DEFAULT_COUNT_DT )
    else:
        result = np.zeros( number_of_references, 
                           dtype = DEFAULT_COUNT_DT )

    for i in range(lower_relative_index, upper_relative_index + 1):
        result += dataset[i*number_of_references:\
                          (i+1)*number_of_references]

    return result

############################################################    

def group_individual_dataset_list(dataset, 
                                  number_of_references,
                                  range_lower,
                                  range_upper,
                                  lower_relative_index,  
                                  upper_relative_index):
    """
    Grabs a slice of the dataset for a given relative index range
    Adds the corresponding read lengths in the first column.
    In contrast to the sum_individual_dataset_list function,
    it does not add count values up.


    Parameters
    ----------

    dataset: Actual dataset of count values from hdf5

    number_of_references: number of references (transcripts)
                          in the annotation of the ribo file

    range_lower: actual minimum read length to be taken
                 This is the initial value of the first column
                 of the resulting np array

    range_upper: actual maximum read length to be taken
                 This is the terminal value of the first column
                 of the resulting np array

    lower_relative_index: relative position of the initial
                          length relative to the minimum 
                          read length of the ribo file
                          i.e.: length_range_min - ribo_length_min

    upper_relative_index: relative position of the terminal
                          length relative to the maximum 
                          read length of the ribo file
                          i.e.: lower_relative_index + \
                                (range_max - range_min)

    Returns
    -------
    A numpy array of the following  form where
    the first column corresponds to the read lengths
    that the counts are coming from
    [ range_lower, count_values_of_length_range_lower ]
    [ range_lower+1, count_values_of_length_range_lower+1 ]
    [ ..........., count_values_.... ]
    [ range_upper, count_values_of_length_range_upper ]

    """
    
    length_column = np.repeat( range(range_lower, range_upper + 1), 
                               number_of_references )
    length_column = np.array( length_column , ndmin = 2, 
                              dtype = DEFAULT_COUNT_DT )

    data_in_range_of_interest = \
        dataset[ lower_relative_index * number_of_references :\
                 (upper_relative_index + 1) * number_of_references ]


    return np.concatenate( (length_column.T, 
                               data_in_range_of_interest), 
                                 axis = 1 )

#######################################################

def aggregate_by_length( ribo_handle, 
                         dataset_name,
                         range_lower, 
                         range_upper,
                         sum_values      = False, 
                         experiment_list = []  ):

    """
    Aggregates values of a dataset by grouping by length 
    or by summing them up for a given read length range.

    This is done for every experiment in the given experiment
    in the experiment_list.
    If the list is empty, then all experiments in ribo
    file handle are taken.

    Parameters
    ----------
    ribo_handle : HDF5 handle of the open ribo file

    dataset_name : Used as a path for the dataset of interest
                   in the ribo file

    range_lower : initial read length to be taken

    range_upper : initial read length to be taken

    sum_values : If true, values are summed up
                 in the given read length range.
                 If false, all the summands are reported
                 with their corresponding read length

    experiment_list : The experiments to be reported.
                  If the list is empty, all experiments in the ribo file are to be taken

    Returns
    -------
    A 2D Numpy Array consisting of Aggregated data by length

    """
    
    if experiment_list != []:
        actual_experiment_list = experiment_list
    else:
        actual_experiment_list = get_experiment_names( ribo_handle )

    dataset_list = tuple(map( lambda x: \
                          ribo_handle[ EXPERIMENTS_name ][x][dataset_name][...], 
                          actual_experiment_list ) )

    read_len_min, read_len_max = get_read_length_range(ribo_handle)

    lower_relative_index = range_lower - read_len_min
    if lower_relative_index < 0:
        raise Exception("Lower range ({}) can not be".format(range_lower) + \
                        " smaller than minimum" +\
                        " read length ({})".format(read_len_min))

    upper_relative_index = (range_upper - range_lower) + lower_relative_index  
    if range_upper < range_lower or range_upper > read_len_max:
        raise Exception("Upper range must be less than the maximum"
                        " read length of the ribo and" 
                        "greater than the lower range")
    
    number_of_references = get_number_of_references(ribo_handle)

    if sum_values:
        aggregate_func = partial(sum_individual_dataset_list,
                          number_of_references = number_of_references, 
                          lower_relative_index = lower_relative_index, 
                          upper_relative_index = upper_relative_index)
    else:
        aggregate_func = partial(group_individual_dataset_list, 
                          number_of_references = number_of_references, 
                          range_lower          = range_lower,
                          range_upper          = range_upper,
                          lower_relative_index = lower_relative_index,  
                          upper_relative_index = upper_relative_index)

    aggregated_data = tuple( map(aggregate_func, dataset_list ) )

    return aggregated_data

############################################################

def aggregate_start_site(ribo_handle, range_lower, range_upper,
                         sum_values      = False,
                         experiment_list = []):
    """
    A wrapper function for aggregate_by_length for start sites
    See the documentation of aggregate_by_length for details
    """

    dataset_name = METAGENE_name + "/" + REF_DG_START_SITE_COV
    return aggregate_by_length(ribo_handle, dataset_name,
                         range_lower, range_upper,
                         sum_values, experiment_list)

def aggregate_stop_site(ribo_handle, range_lower, range_upper, 
                        sum_values      = False,
                        experiment_list = []):
    """
    A wrapper function for aggregate_by_length for stop sites
    See the documentation of aggregate_by_length for details
    """

    dataset_name = METAGENE_name + "/" + REF_DG_STOP_SITE_COV
    return aggregate_by_length(ribo_handle, dataset_name,
                         range_lower, range_upper, 
                         sum_values, experiment_list)

def aggregate_region_counts( ribo_handle, range_lower, range_upper, 
                             sum_values      = False,
                             experiment_list = [] ):
    """
    A wrapper function for aggregate_by_length for region counts
    See the documentation of aggregate_by_length for details
    """

    dataset_name = REF_DG_REGION_COUNTS + "/" + REF_DG_REGION_COUNTS
    return aggregate_by_length(ribo_handle, dataset_name,
                         range_lower, range_upper, 
                         sum_values, experiment_list)
