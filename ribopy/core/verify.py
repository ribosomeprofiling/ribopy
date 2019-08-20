# -*- coding: utf-8 -*-

import functools

import h5py 

from .get_gadgets import *

from .exceptions import ExperimentDoesntExist


#####################################################################
## TODO ?? Revise this function
def prompt_user(message, attr_exists, force):
    if not attr_exists or force:
        return True
    else:
        print(message)
        user_response = input("Do you want to continue? y/n  : ")

        if user_response.lower().startswith("y"):
            return True
        else:
            exit(1)
    


def handle_error(error_message, explanation_message):
    if not error_message:
        return 0
    else:
        error_message += "\n" + explanation_message
        raise ValueError(error_message)

"""
def verify_read_length_range( ribo_handle, 
                              range_lower, range_upper ):
    
    error_message = ""
    length_min, length_max = get_read_length_range(ribo_handle)
    explanation_message    = " ".join(["This ribo file has min read length:",
                                     str(length_min), "and max read length:", 
                                     str(length_max), "."  ])

    if range_lower < length_min or range_lower > length_max:
        error_message = \
          "Given lower range {} is not in the correct interval.".\
             format(range_lower) 

    if range_upper > length_max or range_upper < length_min:
        error_message = \
          "Given upper range {} is not in the correct interval.".\
             format(range_upper)

    return handle_error(error_message, explanation_message) 
"""

"""

def verify_experiment_names (ribo_handle, experiment_names):

    existing_experiment_names = get_experiment_names(ribo_handle)
    error_message = ""
    explanation_message = "Existing experiments are: \n"+\
                          ", ".join(existing_experiment_names)

    for s in experiment_names:
        if s not in existing_experiment_names:
            error_message += \
              "The experiment {} does not exist.\n".format(s)

    return handle_error(error_message, explanation_message)
"""

def verify_extended_region_name(region_name):
    if region_name not in EXTENDED_REGION_names:
        error_message = "No region exists with the name " + region_name
        error_message += " Available regions are:" + \
                          join(", ".EXTENDED_REGION_names)
        raise ValueError(error_message)

def verify_site_type(site_type):
    valid_site_types = (REF_DG_START_SITE_COV, REF_DG_STOP_SITE_COV)
    if site_type not in valid_site_types:
        error_message = "You provided {} as site type.\n".format(site_type)
        error_message+= "Site type can be either {} or {}.".format(
                          *valid_site_types)
        raise ValueError(error_message)


###########################################################
##  D E C O R A T O R S
###########################################################


def make_cli_function(func):
    """
    This function decorator catches RiboPy specific errors and
    prints them on the standard error
    """
    @functools.wraps(func)
    def cli_func_wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except RiboBaseError as e:
            exit(e)
    return cli_func_wrapper


def check_experiment_list_in_ribo_handle(func):
    """
    NOTE THAT THIS DECORATOR IMPOSES CERTAIN
    NAMING CONVENTIONS ON THE ARGUMENTS OF func.

    If there are any experiments that does not exist in the
    ribo, raises an exception.

    Assumptions:
    The first argument of the decorated fucntion is ribo handle
    or it has ribo_handle as a keyword argument

    Experiemnts are given as a list in "experiemnt_list"
    or one experiment as 
    """
    @functools.wraps(func)
    def check_experiments_wrapper(*args, **kwargs):

        # Ribo file is either the first argument
        # or given in the keyword argument entitled "ribo_file"
        ribo_handle = kwargs.get("ribo_handle", None)
        if not ribo_handle:
            ribo_handle = args[0]

        # Either there is an experiment_list
        # or experiment_name
        # or just name in the keyword arguments
        # so we construct experiment list accordingly
        experiment_list       = []
        ribo_experiments      = []
        not_found_experiments = []

        if kwargs.get("experiment_list", None):
            experiment_list = kwargs["experiment_list"]
        elif kwargs.get("experiment_name", None):
            experiment_list = [ kwargs["experiment_name"] ]
        elif kwargs.get("name", None):
            experiment_list = [ kwargs["name"] ]

        if len(experiment_list) > 0:
            ribo_experiments = get_experiment_names(ribo_handle)

        for e in experiment_list:
            if e not in ribo_experiments:
                not_found_experiments.append(e)

        if len(not_found_experiments) > 0:
            plural_1, plural_2 = "", ""

            if len(not_found_experiments) > 1 :
                plural_1 = "s"
                plural_2 = "es" 

            exp_str = ", ".join(not_found_experiments)
            error_message = "The following experiment{} do{} not exist:\n"\
                               .format(plural_1, plural_2) + \
                               exp_str
            raise ExperimentDoesntExist(error_message)

        return func(*args, **kwargs)

    return check_experiments_wrapper


def check_length_range(func):
    """
    Compares the boundaries of the arguments
    to the minimum and maximum length of the ribo file

    If these ranges are not "compatible",
    raises an Exception.

    Assumption:
    -----------
    kwargs of the function decorated must have the keys
    range_lower and range_upper
    """
    @functools.wraps(func)
    def check_length_range_wrapper(*args, **kwargs):
        ribo_handle = kwargs.get("ribo_handle", None)
        if not ribo_handle:
            ribo_handle = args[0]  

        min_read_length, max_read_length = \
            get_read_length_range(ribo_handle)


        # If either of the boundaries are not set,
        # set them to ribo length boundaries.
        if not kwargs["range_lower"] :
            kwargs["range_lower"] = min_read_length
        if not kwargs["range_upper"] :
            kwargs["range_upper"] = max_read_length


        range_lower = kwargs["range_lower"]
        range_upper = kwargs["range_upper"]

        error_message = ""

        # We don't expect this but we check it just in case.
        if max_read_length < min_read_length:
            error_message = "There is a problem with the read length" + \
                            " range of this ribo file. Minimum read length " + \
                            "is greater than maximum length."
            error_message += "Min Read Length = {},\n Max Read Length = "\
                .format(min_read_length, max_read_length)

        if range_lower > range_upper:
            error_message = "Range Lower = {} can not be larger than " + \
                            "Range Upper = {}".format(range_upper, range_upper)

        if range_lower < min_read_length:
            error_message = "Minimum Length can not be less than {}"\
                              .format(min_read_length)

        if range_upper > max_read_length:
            error_message = "Maximumu read length can not be more than {}"\
                              .format(max_read_length)

        if error_message:
            raise InvalidLengthRange(error_message)
        else:
            return func(*args, **kwargs)

    return check_length_range_wrapper
