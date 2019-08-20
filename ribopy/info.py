# -*- coding: utf-8 -*-

import h5py

from .settings import *
from .core.get_gadgets import *
from .io.metadata import get_ribo_file_metadata_for_display,\
                         metadata_dict_to_aligned_str

from .core.verify import make_cli_function
###########################################################

@make_cli_function
def ribo_file_info(ribo_file):
    """
    Returns metadata about the ribo file
    and / or a given experiment if given.

    Parameters
    ----------

    ribo_file : h5py.File
        Path of the ribo file

    """
    with h5py.File(ribo_file, "r") as ribo_handle:
        ribo_info_str = ribo_info(ribo_handle)

    return ribo_info_str



def ribo_info_base(ribo_handle):
    ribo_metadata = get_ribo_file_metadata_for_display(ribo_handle)

    display_lines = metadata_dict_to_aligned_str(ribo_metadata)

    return "Ribo File Info:\n" +\
           "---------------\n" + display_lines



def _format_experiment_info( exp_info_list ):
    """
    Formats experiment information by aligning columns

    The list exp_info_list is of the form
    [ [name, reads, rnaseq, metadata] ]
    """
    num_of_columns = 5
    col_max_lens   = []
    spacer         = "  "

    for i in range(num_of_columns):
        this_max = 0
        for entry_list in exp_info_list:
            this_length = len(entry_list[i])
            if this_length > this_max:
                this_max = this_length
        col_max_lens.append(this_max)

    formatted_list = []
    for entry_list in exp_info_list:
        raw_name  = entry_list[0]
        line_list = [ raw_name + \
                      (" " * (col_max_lens[0] - len(raw_name))) ]

        for i in range(1,num_of_columns):
            raw_val = entry_list[i]
            val     = ( " " * (col_max_lens[i] - len(raw_val)) ) + \
                        raw_val

            line_list.append(val)

        formatted_list.append(spacer.join( line_list ))

    return formatted_list


def ribo_info(ribo_handle):
    """
    Returns essential information about the ribo file.
    """

    experiment_list = get_experiment_names(ribo_handle)
    display_lines    = []
    display_base     = ribo_info_base(ribo_handle)

    experiment_str = "Experiments:\n------------" \
        if len(experiment_list) > 1 else "Library:\n--------"
    display_lines.append("\n" + experiment_str + " ")

    raw_exp_info_list = [[ "Name", "Reads", "Coverage", "RNA-Seq", "Metadata" ]]

    for experiment in experiment_list:
        raw_exp_info_list.append( get_experiment_info(ribo_handle, experiment) )


    display_lines += _format_experiment_info(raw_exp_info_list)

    return display_base + "\n" + "\n".join(display_lines)


def get_experiment_info(ribo_handle, experiment):
    """
    Returns essential information about the experiment.

    The returned value is a list of the form
    [name, reads, coverage, rnaseq, metadata]
    """
    exp_handle = ribo_handle[EXPERIMENTS_name][experiment]

    if has_coverage_data(ribo_handle, experiment):
        coverage_str = "*"
    else:
        coverage_str = ""

    if rnaseq_exists(ribo_handle, experiment):
        rnaseq_str = "*"
    else:
        rnaseq_str = ""

    if exp_handle.attrs.get(USER_METADATA, ""):
        metadata = "*"
    else:
        metadata = ""

    total_reads = get_total_reads(ribo_handle, experiment)

    return [experiment, str(total_reads), coverage_str, rnaseq_str, metadata]
