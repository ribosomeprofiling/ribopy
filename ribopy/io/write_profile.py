# -*- coding: utf-8 -*-

from collections import OrderedDict
import numpy as np
import pandas as pd

from ..settings import *

###########################################################

def _merge_profile( profile_list ):
    array_size = len(profile_list[0])

    for a in profile_list[1:]:
        assert array_size == len(a)

    return np.concatenate( profile_list ) 

def _merge_region( region_counts, 
                   length_min, length_max ):
    pass


def _write_metagene( metagene_handle, 
                     start_site, stop_site, 
                     metagene_radius ):

    metagene_handle.create_dataset(
                            REF_DG_START_SITE_COV,
                            data        = start_site,
                            dtype       = SITE_COVERAGE_DT,
                            compression = DEFAULT_COMPRESSION, 
                            fletcher32  = DEFAULT_FLETCHER32)

    metagene_handle.create_dataset(
                            REF_DG_STOP_SITE_COV,
                            data        = stop_site,
                            dtype       = SITE_COVERAGE_DT,
                            compression = DEFAULT_COMPRESSION, 
                            fletcher32  = DEFAULT_FLETCHER32)

    metagene_handle.attrs[METAGENE_RADIUS_name] = metagene_radius


def _write_region_counts( region_handle, 
                          region_counts, 
                          left_span, right_span ):

    region_handle.create_dataset(
                            REF_DG_REGION_COUNTS,
                            data        = region_counts,
                            dtype       = DEFAULT_COUNT_DT,
                            compression = DEFAULT_COMPRESSION, 
                            fletcher32  = DEFAULT_FLETCHER32)

    region_handle.attrs[LEFT_SPAN_name]  = left_span
    region_handle.attrs[RIGHT_SPAN_name] = right_span


#####################################################################

def _write_coverage(exp_handle, coverages, store_coverage):
    """
    ??? Incomplete documentation
    """

    # First coverage is for the first length
    first_coverage    = np.concatenate( tuple( coverages[0].values() ) )
    # The length of first coverage gives us the total
    # number of nucleotides in this transcript set
    total_tcript_nucs = len(first_coverage)

    total_coverage_size = len(coverages) * len(first_coverage)

    # allocating space for the whole array first and thenm filling it up
    # should be more efficient than concatanating arrays one after the other
    # though we might need to optimize this part  
    coverage_np         = np.zeros( total_coverage_size, 
                                    dtype = TRANSCRIPT_COVERAGE_DT )

    coverage_np[0:total_tcript_nucs] = first_coverage

    x = total_tcript_nucs
    y = x + total_tcript_nucs

    for c in coverages[1:]:
        this_coverage     = np.concatenate( tuple( c.values() ) )
        coverage_np[x:y]  = this_coverage
        x  = y
        y += total_tcript_nucs 

    if store_coverage:
        exp_handle.create_group(REF_DG_COVERAGE)
        coverage_handle = exp_handle[REF_DG_COVERAGE] 

        coverage_handle.create_dataset(
                            REF_DG_COVERAGE,
                            data        = coverage_np, 
                            dtype       = TRANSCRIPT_COVERAGE_DT,
                            compression = DEFAULT_COMPRESSION, 
                            fletcher32  = DEFAULT_FLETCHER32)

    return coverage_np.sum()



#####################################################################

def write_profile( exp_handle , profiles,
                   length_min , length_max,
                   left_span,   right_span, 
                   metagene_radius,
                   store_coverage = False ):
    
    start_site_profiles = list( map( lambda x: x[REF_DG_START_SITE_COV], 
                                     profiles ) )
    stop_site_profiles  = list( map( lambda x: x[REF_DG_STOP_SITE_COV], 
                                     profiles ) )
    region_counts       = list( map( lambda x: x[REF_DG_REGION_COUNTS], 
                                     profiles ) )
    coverages           = list( map( lambda x: x[REF_DG_COVERAGE],
                                     profiles ))

    merged_start_site_profiles = _merge_profile( start_site_profiles )
    merged_stop_site_profiles  = _merge_profile( stop_site_profiles )
    merged_region_counts       = _merge_profile(region_counts)

    exp_handle.create_group(METAGENE_name)
    metagene_handle = exp_handle[METAGENE_name]
    _write_metagene( metagene_handle, 
                     merged_start_site_profiles, 
                     merged_stop_site_profiles,
                     metagene_radius )

    exp_handle.create_group(REF_DG_REGION_COUNTS)

    region_handle = exp_handle[REF_DG_REGION_COUNTS]

    _write_region_counts( region_handle, 
                          merged_region_counts,
                          left_span, right_span )

    total_reads = _write_coverage(exp_handle, coverages, 
                                  store_coverage)

    exp_handle.attrs[ATTRS_TOTAL_READS] = total_reads

    
    

        
    
