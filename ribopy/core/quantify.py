# -*- coding: utf-8 -*-

import numpy as np
from .coverage import find_coverage
from .metagene import find_site_coverage
from .region_counts import find_region_counts
from ..settings import *
from ..io.file import flex_open 


################################################################

def quantify_experiment( 
        input_reads_file, 
        ref_names, 
        ref_lengths, 
        region_coordinates, 
        metagene_radius, 
        left_span, 
        right_span):
    """
    Does all essential quantification for an input bed file.
    The results are packed in a dictionary called "essential_profile".
    The input reads can come from any bed file. 
    Though in riboflow, when we call this function,
    the input reads are of a fixed length.

    Parameters
    ----------
    input_reads : str 
        In bed format. It can also be an open stream.

    ref_names   : list(str) 
        Reference names in an array.

    ref_lengths : list(int)
        Reference lengths with the same order as names.

    region_coordinates : list( [start_u5, stop_u5], 
                               [start_cd, stop_cd]
                               [start_u3, stop_u3])
                         
                         Annotation of the reference.
                         In triplets of the form [start, stop]
                         for UTR5, CDS and UTR3
    metagene_radius : int
                      nucleotides to the left / right of start / stop site
                      for metagene analysis

    left_span: int
               Nucleotides to the left of start / stop site for region counts

    left_span: int
               Nucleotides to the right of start / stop site for region counts


    Returns
    -------
    essential_profile : An array of 

    """

    with flex_open( input_reads_file, "r" ) as input_stream:
        coverage = find_coverage(  input_stream , ref_names, ref_lengths )

    essential_profile = dict()

    essential_profile[REF_DG_COVERAGE] = coverage    

    essential_profile[REF_DG_START_SITE_COV] = \
        find_site_coverage(coverage, metagene_radius, 
             region_coordinates, site_type = "start")

    essential_profile[REF_DG_STOP_SITE_COV] = \
        find_site_coverage(coverage, metagene_radius, 
             region_coordinates, site_type = "stop")

    essential_profile[REF_DG_REGION_COUNTS] = \
        find_region_counts( coverage, 
             region_coordinates, left_span, right_span )

    return essential_profile






