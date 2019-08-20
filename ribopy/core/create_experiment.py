# -*- coding: utf-8 -*-

from functools import partial
from multiprocessing import Pool
import os

from ..io.separate_by_length import separate_by_length
from .quantify import quantify_experiment
from .get_gadgets import *
from ..io.write_profile import write_profile
from ..metadata import set_metadata
from ..settings import *


def create_experiment(ribo_exp_handle,
                   experiment_name, 
                   alignment_file_handle,
                   ref_names, 
                   ref_lengths, 
                   region_coordinates,
                   metagene_radius, 
                   left_span,  right_span,
                   length_min, length_max,
                   metadata        = None, 
                   store_coverage  = False,
                   nprocess        = 1,
                   tmp_file_prefix = "" ):

    """
    ??? Missing documentation ???
    """   

    # Step 1: Write Metadata
    if metadata:
        set_metadata(ribo_exp_handle, "", metadata = metadata)

    

    # Step 2: Separate files by length
    if tmp_file_prefix == "":
        tmp_file_prefix = experiment_name + "."

    files_by_length = separate_by_length(
                          input_stream = alignment_file_handle, 
                          length_max   = length_max, 
                          length_min   = length_min, 
                          file_prefix  = tmp_file_prefix)

    # files_by_length is a tuple of the form
    #   ( (file_1.bed, length_min), (file_2.bed, length_min+1), ... , 
    #     (file_k.bed, length_max) )


    quantify = partial( quantify_experiment, 
                           ref_names          = ref_names, 
                           ref_lengths        = ref_lengths, 
                           region_coordinates = region_coordinates, 
                           metagene_radius    = metagene_radius, 
                           left_span          = left_span, 
                           right_span         = right_span )

    # Step 3: Map these files to count structures ( metagene, regions counts etc )
    #quant_data_by_length = map( quantify, files_by_length )
    with Pool( nprocess ) as p:
        profiles = tuple( p.map(quantify, files_by_length) )

    for f in files_by_length:
        os.remove(f)

    # Step 4: write the result to the ribo file
    # The results are merged in the write function
    write_profile( exp_handle      = ribo_exp_handle , 
                   profiles        = profiles,
                   length_min      = length_min, 
                   length_max      = length_max,
                   left_span       = left_span,   
                   right_span      = right_span, 
                   metagene_radius = metagene_radius,
                   store_coverage  = store_coverage )
    