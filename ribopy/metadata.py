# -*- coding: utf-8 -*-
from collections import OrderedDict
import os
from io import IOBase

import h5py
import yaml

from .settings import *
from .io.file import get_alignment_file_handle
from .core.get_gadgets import *
from .core.exceptions import *
from .io.metadata import metadata_dict_to_aligned_str
from .io.file import (open_by_extension,
                      read_all_lines, 
                      flex_open) 
from ._version import (__format_version__,
                      __version__)
from .core.verify import prompt_user, make_cli_function ,\
                         check_experiment_list_in_ribo_handle


###########################################################

# From : https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
# This is a fix to ordering metadata
# Though it is not very important, we want to see the
# metadata in the same order.

_mapping_tag = yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG

def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())

def dict_constructor(loader, node):
    return OrderedDict(loader.construct_pairs(node))

yaml.add_representer(OrderedDict, dict_representer)
yaml.add_constructor(_mapping_tag, dict_constructor)

# TODO
# Check experiment existence

###########################################################

@check_experiment_list_in_ribo_handle
def set_metadata(ribo_handle, name, metadata):
    """
    Set the user-provided metadata

    Note that this will over-write the existing metadata

    Metadata needs to be in yaml format.
    It can be a dictionary, string or a filestream.

    If no name is given, the metadata of the root ribofile
    is going to be set.

    Parameters
    ----------
    ribo_handle: h5py.File
        Open ribo file handle

    name: str
        Name of the experiment

    metadata: dict, str or IOStream
        metadata in yaml format

    Returns
    -------
    error: str
      Error returned by yaml module.
      It is the empty string if there are no errors.
    """

    handle = ribo_handle

    if name:
        handle = ribo_handle[EXPERIMENTS_name][name]

    metadata_type = type(metadata)
    # make sure the metadata is fine
    # before storing
    if metadata_type == str or issubclass(metadata_type, IOBase):
        try:
            meta_dict = yaml.safe_load(metadata)
            if not issubclass(type(meta_dict), dict):
                e_msg = " Entries must be in the form\n" +\
                        " key : value"
                raise (yaml.error.YAMLError(e_msg))
        except yaml.error.YAMLError as e:
            return "Invalid Metadata\n" + str(e)
    else:
        meta_dict = metadata

    meta_str = yaml.safe_dump(meta_dict)

    handle.attrs[USER_METADATA] = meta_str

    return ""
 
@make_cli_function
def set_metadata_wrapper(ribo_file, name, meta_file, force):
    """
    Wrapper for set_metadata
    """

    with flex_open(meta_file) as file_handle,\
          h5py.File(ribo_file, "r+") as ribo_handle:

        # <Prompt>
        handle = ribo_handle

        if name:
            handle = ribo_handle[EXPERIMENTS_name][name]

        prompt_attr    = handle.attrs.get(USER_METADATA, None)
        prompt_message = "This will over-write the existing metadata."
        prompt_user(message     = prompt_message, 
                    attr_exists = prompt_attr , 
                    force       = force)

        # </Prompt>


        error = set_metadata(ribo_handle = ribo_handle,
                             name        = name,
                             metadata    = file_handle)

        if error:
            print(" Failure:\n", error)

@check_experiment_list_in_ribo_handle
def get_metadata(ribo_handle, name):
    """
    Returns user-defined metadata in a dictionary

    If no name is given, the metadata of the root ribo file
    is returned.

    Parameters
    ----------
    ribo_handle: h5py.File
        Open ribo file handle

    name: str
        Name of the experiment

    Returns
    -------
    metadata in dictionary form    
    """

    handle = ribo_handle

    if name:
        handle = ribo_handle[EXPERIMENTS_name][name]

    metadata = handle.attrs.get(USER_METADATA, "")

    return yaml.safe_load(metadata)

@make_cli_function
def get_metadata_wrapper(ribo_file, name):
    """
    Wrapper for get_metadata
    """

    with h5py.File(ribo_file, "r") as ribo_handle:
        metadata = get_metadata(ribo_handle, name)
    
    if metadata: 
        print(yaml.safe_dump(metadata, default_flow_style = False))

    # We can not align nested dictionaries for now.
    # So we are not using fancy output
    #print( metadata_dict_to_aligned_str(metadata) )

@check_experiment_list_in_ribo_handle
def delete_metadata(ribo_handle, name):
    """
    Deletes user-provided metadata

    Parameters
    ----------
    ribo_handle: h5py.File
        Open ribo file handle

    name: str
        Name of the experiment

    """

    handle = ribo_handle

    if name:
        handle = ribo_handle[EXPERIMENTS_name][name]

    if handle.attrs.get(USER_METADATA, None):
        print("Deleting metadata")
        del handle.attrs[USER_METADATA]

    return ""

@make_cli_function
def delete_metadata_wrapper(ribo_file, name, force):
    """
    Wrapper for delete_metadata
    """

    with h5py.File(ribo_file, "r+") as ribo_handle:

        # <Prompt>
        handle = ribo_handle

        if name:
            handle = ribo_handle[EXPERIMENTS_name][name]

        prompt_attr    = handle.attrs.get(USER_METADATA, None)
        prompt_message = "Existing metadata is going to be deleted."
        prompt_user(message     = prompt_message, 
                    attr_exists = prompt_attr , 
                    force       = force)

        # </Prompt>
        delete_metadata(ribo_handle, name)
        
@check_experiment_list_in_ribo_handle
def has_metadata(ribo_handle, name = None):
    """
    Checks if the ribo file or the experiment has metadata
    """
    if name:
        this_handle = ribo_handle[EXPERIMENTS_name][name]
    else:
        this_handle = ribo_handle
        
    if this_handle.attrs.get(USER_METADATA, None):
        return True
    else:
        return False
    
