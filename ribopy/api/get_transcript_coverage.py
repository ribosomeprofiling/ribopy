# -*- coding: utf-8 -*-

from functools import partial, wraps
from collections import OrderedDict

import h5py
import pandas as pd
import numpy as np

from ..core.exceptions import *
from ..settings import *

def check_index(ribo, range_lower, range_upper):
    if range_upper - range_lower < 0:
        raise InvalidLengthRange(
         "range_lower = {} must be less than or equal to range upper {}"\
         .format( range_lower, range_upper ))
    
    if range_lower < ribo.minimum_length :
        raise InvalidLengthRange(
           "range_lower can not be less than the minimum length of the ribo (={})"\
               .format(ribo.minimum_length))
               
    if range_upper > ribo.maximum_length :
        raise InvalidLengthRange(
           "range_upper can not be more than the maximum length of the ribo (={})"\
               .format(ribo.maximum_length))

def _get_transcript_coverage(ribo,
                            experiment_name,
                            transcript,
                            transcript_length,
                            total_length,
                            range_lower = 0, 
                            range_upper = 0,
                            sum_lengths = False):
    """
    Get the coverage vector of a single transcript.    
    """
    
    coverage_handle = ribo._handle[EXPERIMENTS_name][experiment_name]\
                         [REF_DG_COVERAGE][REF_DG_COVERAGE]
    
    if range_upper == 0:
        range_upper = ribo.maximum_length
    if range_lower == 0:
        range_lower = ribo.minimum_length
    
    check_index(ribo, range_lower, range_upper)
    
    range_size      = (range_upper - range_lower) + 1
    coverage_series = np.zeros(shape = ( range_size, transcript_length ), 
                               dtype = TRANSCRIPT_COVERAGE_DT)
                               
    main_offset = (range_lower - ribo.minimum_length) * total_length
    transcript_offset = ribo.transcript_offsets[transcript]                               
    
    for i in range(range_lower, range_upper + 1):
        array_index                  = i - range_lower
        slice_start                  = main_offset + transcript_offset
        slice_end                    = slice_start + transcript_length
        coverage_series[array_index] = coverage_handle[slice_start : slice_end]
        main_offset                 += total_length
    
    coverage_df = pd.DataFrame(coverage_series, 
                               index = tuple( range(range_lower, range_upper + 1)),
                                )
    coverage_df.index.name = "length"
    
    if sum_lengths:
        coverage_df = pd.DataFrame( coverage_df.sum(axis = 0) ).T

    return coverage_df
