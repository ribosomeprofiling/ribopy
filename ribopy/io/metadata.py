# -*- coding: utf-8 -*-

from collections import OrderedDict
import datetime

import yaml

from ..settings import *
from ..core.get_gadgets import get_experiment_names
from .file import flex_open

##################################################3


def metadata_dict_to_aligned_str(metadata_dict):
    """
    Aligns metadata keys.

    Returns the resulting string.
    """

    if not metadata_dict or len(tuple( metadata_dict.keys() )) == 0:
        return ""

    max_key_len = max( map(len, metadata_dict.keys()) )
    
    def aligner(x):
        return  x + ( " " * ( max_key_len - len(x) ) )\
                  + " : " + str(metadata_dict[x])
    
    aligned_keys = tuple(map( aligner, metadata_dict.keys()))
    
    return "\n".join(aligned_keys)


def read_metadata_from_yaml(metadata):
    with flex_open(metadata) as metadata_input:
        result_dict = yaml.load(metadata_input, Loader=yaml.FullLoader)
    return result_dict
   

def get_ribo_file_metadata_for_display(ribo_handle):
    ribo_metadata = OrderedDict()

    for key, value in RIBO_METADATA_FOR_DISPLAY.items():
        ribo_metadata[value] = ribo_handle.attrs[key]

    # format time
    time_in_seconds = ribo_metadata[RIBO_METADATA_FOR_DISPLAY[ATTRS_TIME]] 
    val = datetime.datetime.fromtimestamp(time_in_seconds)
    ribo_metadata[RIBO_METADATA_FOR_DISPLAY[ATTRS_TIME]] = \
        val.strftime('%Y-%m-%d %H:%M:%S') 

    return ribo_metadata
