# -*- coding: utf-8 -*-

import h5py
import pandas as pd
import numpy as np

from ..settings import *
from .exceptions import *

###########################################################

def get_reference_names(h5_handle, alias = None):
    """
    Returns reference (transcript) names as a string array 
    
    If "alias" is provided, it returns the alias instead.
    """
    original_names =  h5_handle[REFERENCE_name]\
                          [REF_DG_REFERENCE_NAMES][...].astype(str)
    
    if alias is None:
        return original_names
    else:
        alias_length = len(alias)
        alias_type   = type(alias)
        if alias_length != len(original_names):
            raise AliasError("alias length must be {}"\
                            .format(len(original_names)))
        if alias_type != list and alias_type != tuple:
            raise AliasError("Alias my be an array of tuple of strings!")
        return alias



def get_reference_lengths(h5_handle):
    """
    Returns reference lengths in the same order as
    reference names (get_reference_names)
    """
    return h5_handle[REFERENCE_name]\
                    [REF_DG_REFERENCE_LENGTHS][...].astype(int)

def get_number_of_references(h5_handle):
    """
    Returns the total number of transcripts (references)
    Use to create the ribo file
    """
    return len(h5_handle[REFERENCE_name]\
                        [REF_DG_REFERENCE_NAMES][...])

def get_region_boundaries(h5_handle):
    """
    Returns an array of pairs where pairs are of the form
    ( [UTR5_start=0, UTR5_end], [CDS_start, CDS_end], [UTR3_start, UTR3_end] )
    The coordinates are in bed format.
    Therefore it is zero based, first value is inclusive 
    and the second value is exclusive.
    """
    annotation_matrix = h5_handle[REFERENCE_name][REF_ANNOTATION_NAME][...]
    boundaries = map( lambda x: [(0, x[0]) , (x[0], x[1]), (x[1], x[2])] , 
                      annotation_matrix )

    return tuple(boundaries)

def get_experiment_names(ribo):
    """
    Returns experiment names in a ribo file.
    Initially, a ribo file has only one experiment
    But ribo files can be merged to contain more than one experiment.
    """
    experiments_handle = ribo[EXPERIMENTS_name]
    return tuple( experiments_handle.keys() )


def get_read_length_range(ribo):
    """
    Returns the minimum and maximum read lengths
    """
    return ( int(ribo.attrs[LENGTH_MIN_name]), 
             int(ribo.attrs[LENGTH_MAX_name]) )

def get_metagene_radius(ribo):
    """
    Returns the number of nucleotides to the left and right of
    start / stop sites in metagene analysis tables.
    """
    return ribo.attrs[ATTRS_METAGENE_RADIUS]

def get_left_right_spans(ribo):
    """
    In counting number of reads going to each region (UTR5, CDS and UTR3),
    We artificially introduce two regions

    UTR5_junction: A region between UTR5 and CDS containing the start site
                   [start_site - left_span, start_site + leftspan]

    UTR3_junction: A region between CDS and UTR3 containing the stop site
                   [stop_site - left_span, stop_site + leftspan]

    This function returns a tuple of the form:
    (left_span, right_span)

    """
    return (ribo.attrs[ATTRS_LEFT_SPAN],\
            ribo.attrs[ATTRS_RIGHT_SPAN])


def has_coverage_data(ribo, experiment):
    """
    Tells if the given experiment has coverage data
    stored.
    Note that storage of coverage data is optional so
    it may not exist in all experiments.
    """

    exp_handle = ribo[EXPERIMENTS_name][experiment]

    # Note that coverage exists  inside a datagroup
    # So we check the existince of the datagroup 
    # and the dataset 
    if REF_DG_COVERAGE not in exp_handle.keys() or \
        REF_DG_COVERAGE not in exp_handle[REF_DG_COVERAGE].keys():
        return False
    else:
        return True

def get_total_reads(ribo, experiment):

    exp_handle = ribo[EXPERIMENTS_name][experiment]

    return exp_handle.attrs[ATTRS_TOTAL_READS]

def rnaseq_exists(ribo, experiment):
    exp_handle = ribo[EXPERIMENTS_name][experiment]

    return  RNASEQ_name in exp_handle.get(RNASEQ_name, [])
