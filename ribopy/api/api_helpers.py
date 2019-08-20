# -*- coding: utf-8 -*-
import functools
from functools import partial

from ..settings import *
from ..io.metadata import get_ribo_file_metadata_for_display
from ..metadata import get_metadata as main_get_metadata
from ..info import ribo_info, get_experiment_info

from ..core.get_gadgets import rnaseq_exists, \
                              has_coverage_data, \
                              get_experiment_names, \
                              has_coverage_data, \
                              rnaseq_exists, \
                              get_total_reads, \
                              get_reference_names

from ..core.exceptions import ExperimentDoesntExist
##################################################################

def _make_ribo_info_dict(ribo_handle):
    """
    Helper function to create a dictionary containing essential
    ribo file attribute.
    """
    
    info_dict                   = get_ribo_file_metadata_for_display(ribo_handle)
    info_dict[EXPERIMENTS_name] = dict() 
    
    experiments = get_experiment_names(ribo_handle)
    for e in experiments:
        # The oupput of get_experiment_infois of the form 
        # [ "Name", "Reads", "Coverage", "RNA-Seq", "Metadata" ]
        exp_info = get_experiment_info(ribo_handle, e)

        info_dict[EXPERIMENTS_name][e] = \
           { "Reads"    : get_total_reads(ribo_handle, e), 
             "Coverage" : has_coverage_data(ribo_handle, e),
             "RNA-Seq"  : rnaseq_exists(ribo_handle, e), 
             "Metadata" : main_get_metadata(ribo_handle, e) }
                                
    return info_dict
    
def _check_experiments(func):
    """
    Checks if the experiments provided in the parameters
    exist in the experiment list of the ribo object.
    
    This decorator assumes that the experiment parameter
    is provided either as keyword parameter or it is the first parameter
    after "self" in the class function definition
    """
    @functools.wraps(func)
    def exp_verifier(*args, **kwargs):
        if kwargs.get("experiments", None):
            experiments = kwargs.get("experiments")
        elif kwargs.get("experiment", None):
            experiments = kwargs.get("experiment")
        elif len(args) > 1:
            experiments = args[1]
        else:
            experiments = None
            
        if not experiments:
            return func(*args, **kwargs)
            
        if type(experiments) == str:
            experiments = [experiments]
        
        valid_experiment_list = args[0].experiments
        for e in experiments:
            if e not in valid_experiment_list:
                error_message = "The experiment '{}' doesn't exist!".format(e)
                raise ExperimentDoesntExist(error_message)
        return func(*args, **kwargs)
    return exp_verifier
