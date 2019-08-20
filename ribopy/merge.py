# -*- coding: utf-8 -*-
from collections import OrderedDict
import time

import h5py
import pandas as pd
import numpy as np
import yaml

from .settings import *
from .core.get_gadgets import *
from .core.metagene import *
from ._version import __format_version__, __version__

##########################################################




#####################################################################
#### C O M P A T I B I L I T Y    F U N C T I O N S
#####################################################################

def check_referenence_compatibility( first_pair, second_pair ):
    """
    ?? Incomplete documentation
    """

    reference_error = "Reference Error\n"

    first_reference_names  = get_reference_names( first_pair[0] )
    second_reference_names = get_reference_names( second_pair[0] )

    if len(first_reference_names) != len(second_reference_names):
        reference_error += "They have different number of references"
        reference_error += "{} != {}".format( len(first_reference_names),
                                              len(second_reference_names))
        return reference_error

    ref_name_comparison = first_reference_names == second_reference_names
    if not ref_name_comparison.all():
        return "Different reference (transcript) names."

    first_reference_lengths  = get_reference_lengths( first_pair[0] )
    second_reference_lengths = get_reference_lengths( second_pair[0] )

    ref_len_comparison = first_reference_lengths == second_reference_lengths
    if not ref_len_comparison.all():
        return "Different reference (transcript) lengths."

    return ""


def check_attribute_compatibility( first_pair, second_pair ):
    """
    ?? Incomplete documentation
    """

    attribute_error = "Attribute Error:\n"

    for attribute in ATTRS_ESSENTIAL_COMPATIBILITY:
        first_attr  = first_pair[0].attrs[attribute]
        second_attr = second_pair[0].attrs[attribute]

        if first_attr != second_attr:
            attribute_error += \
              " The {} atrribute is different:\n".format(attribute)
            attribute_error += "{first} != {second}".format(
                               first  = first_attr,
                               second = second_attr )
            return attribute_error

    return ""


############################################################################

def check_ribo_compatibility_pair( first_pair, second_pair ):
    """
    ?? Incomplete documentation
    """

    first_handle     = first_pair[0]
    first_identifier = first_pair[1]

    second_handle     = second_pair[0]
    second_identifier = second_pair[1]

    error_message = "The ribo files {}, {} are not compatible.\n".format(
                      first_identifier, second_identifier)

    # check attribute compatibility first
    attribute_error = \
        check_attribute_compatibility( first_pair, second_pair )
    if attribute_error:
        raise ValueError( error_message + attribute_error )

    ref_error = check_referenence_compatibility( first_pair, second_pair )
    if ref_error:
        raise ValueError( error_message + ref_error )


def check_if_common_libs_exist( ribo_handle_list ):
    """
    ?? Incomplete documentation
    """

    for i, ribo_1 in enumerate(ribo_handle_list):
        for ribo_2 in ribo_handle_list[i+1:]:
            experiments_1 = set(get_experiment_names(ribo_1[0]) )
            experiments_2 = set(get_experiment_names(ribo_2[0]) )

            common_experiments = experiments_1.intersection(experiments_2)

            if common_experiments :
                identifier_1  = ribo_1[1]
                identifier_2  = ribo_2[1]
                error_message = "The ribos {first} and {second} ".format(\
                                  first = identifier_1,
                                  second = identifier_2)
                error_message += " have common experiments:\n"
                error_message += str(common_experiments)
                raise ValueError(error_message)


def check_ribo_compatibility(ribo_handle_list):
    # Ribo_handle_list is a list of pairs where
    # an element of the list is of the form
    # (hdf5_file_handle, identifier)

    check_if_common_libs_exist(ribo_handle_list)

    for h in ribo_handle_list[1:]:
        check_ribo_compatibility_pair( ribo_handle_list[0], h )


#####################################################################
#### M A I N    F U N C T I O N S
#####################################################################

def initialize( h5_destination_handle, h5_source_handle ):
    """
    ?? Incomplete documentation
    """

    h5_source_handle.copy( REFERENCE_name,
               h5_destination_handle )
    h5_destination_handle.create_group(EXPERIMENTS_name)


###########################################################

def _copy_attributes(h5_destination_handle, h5_source_handle):
    """
    ?? Incomplete documentation
    """

    for key in ATTRIBUTES_TO_BE_COPIED_FOR_MERGE:
        h5_destination_handle.attrs[key] = h5_source_handle.attrs[key]
        
def _copy_ribo_metadata(destination_handle, source_list):
    """
    Tries to merge the metadata of the source_likst into one dictionary
    and write it to destination.
    If source has no metadata, destrination metadata attribute will be set
    but it will be empty.
    """
    merged_metadata_dict = {}
    source_with_metadata = []
    
    for h in source_list:
        if h.attrs.get(USER_METADATA, None):
            metadata_dict = yaml.safe_load(h.attrs[USER_METADATA])
            # Metadata might exists as an emty string
            # So let's make sure it translates to a non-emptry dict.
            if metadata_dict:
                 merged_metadata_dict.update( metadata_dict )
            
    destination_handle.attrs[USER_METADATA] = \
        yaml.safe_dump( merged_metadata_dict )
            
        

def merge_ribos(destination_handle, source_list):
    """
    ?? Incomplete documentation
    """

    # If the sdource list is not coming from pairs,
    # then make it into pairs
    # this way, identifying incompatible ribo files or
    # handles is going to be easier
    if type(source_list[0]) not in (tuple, list):
        source_list = [ ( source_list[i], str(i) ) \
                         for i in range(len(source_list)) ]

    check_ribo_compatibility(source_list)

    source_handle_list = list( map(lambda x: x[0], source_list) )

    initialize(destination_handle, source_handle_list[0])
    _copy_attributes(destination_handle, source_handle_list[0])
    _copy_ribo_metadata(destination_handle, source_handle_list)

    for ribo_handle in source_handle_list:

        for lib in get_experiment_names(ribo_handle):
            exp_path = EXPERIMENTS_name + "/" + lib
            #destination_handle.create_group(exp_path)
            ribo_handle.copy( exp_path, destination_handle[EXPERIMENTS_name] )
            
    destination_handle.attrs[ATTRS_TIME]  = time.time()



def merge_ribo_files(destination_file, source_file_list):
    """
    Merges the experiments in the given source files and writes
    the result in the destination file.

    The ribo files need to be compatible
    (same left / right span, same metagene radius, same reference)

    Because of the compatibility, parameters (attributes), refrence
    etc. of the new file is the same as the merged file,

    The source files are not allowed to have experiments of the same name
    as this creates ambiguity.

    Parameters
    ----------
    destination_file : Destination ribo file path

    source_file_list : List of ribo file paths to be merged

    """

    if len(source_file_list) < 2:
        print("Please provide at least two input ribo files")
        exit(1)

    source_handle_list      = [ (h5py.File(f , "r"), f )
                                for f in source_file_list ]
    destination_ribo_handle = h5py.File( destination_file, "w" )

    merge_ribos( destination_ribo_handle, source_handle_list )

    [s[0].close() for s in source_handle_list]
    destination_ribo_handle.close()
