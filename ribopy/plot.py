# -*- coding: utf-8 -*-

from functools   import partial
from collections import OrderedDict
from functools   import partial, wraps

import h5py
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
# The "use" method is necessary for not seeing 
# gtk error message on the terminal
matplotlib.use('Agg')

from .settings                 import *
from .core.get_gadgets         import *
from .core.aggregate_by_length import aggregate_region_counts,\
                                      aggregate_start_site,\
                                      aggregate_stop_site

from .dump        import *
from .core.verify import verify_extended_region_name,\
                         verify_site_type

from .core.exceptions import MissingExperiment

from .core.verify import make_cli_function,\
                         check_experiment_list_in_ribo_handle,\
                         check_length_range

#####################################################################

def save_plot(plot_func):
    @wraps(plot_func)
    def out_plot_func(*args, **kwargs):
        fig = plt.figure( )
        res = plot_func(*args, **kwargs)
        if "output_file" in kwargs.keys() and kwargs["output_file"]:
            plt.savefig(kwargs["output_file"], bbox_inches = "tight")
        plt.close(fig)
        return fig
    return out_plot_func

def fix_experiment_list(func):
    @wraps(func)
    def exp_fixed_func(*args, **kwargs):
        if not kwargs["experiment_list"]:
            raise MissingExperiment("At least one experiment is required!")
        if type(kwargs["experiment_list"]) == str:
            kwargs["experiment_list"] = [kwargs["experiment_list"]]
        return func(*args, **kwargs)
    return exp_fixed_func


####################################################################    

@check_length_range
@fix_experiment_list
@save_plot
@check_experiment_list_in_ribo_handle
def plot_metagene( ribo_handle, 
                   site_type,
                   *, 
                   experiment_list, 
                   title,
                   range_lower = 0, 
                   range_upper = 0,
                   normalize   = False,
                   output_file = "",
                   colors      = PLOT_COLORS ):
    
    """
    Generates metagene plots at start or stop sites.
    If a file name / path is given, it saves the output to the
    provided file with the given name / path.

    At least one experiment must be provided.
    Unless you provide a custom list of colors,
    at most 7 experiments can be drawn on the same plot
    as there are 7 colors available by default.

    normalize = True divides the coverage vector by the sum of the
    vector (not the whole transcript coverage) and multiplies
    by 100 and gives percentage.
    This is especially useful when plotting experiments together.

    All read lengths in the interval range_lower
    and range_upper are aggregated for the resulting plot.



    Parameters
    ----------
    ribo_handle: h5py.File
       ribo file handle

    site_type: str
       for start site, use "start_site_coverage"
       for stop site use "stop_site_coverage"

    experiment_list: list 
       list of experiment names

    title: str
       Title of the graph

    range_lower: int
       Minimum read length to be taken
       If 0, minimum length of the ribo file is taken

    range_upper: int
       Maximum read length to be taken
       If 0, maximum length of the ribo file is taken

    normalize: Boolean
       If false, plot raw counts
       If True, divide the coverage vector 

    output_file: str
       The path of the output file.
       If it is an NOT empty string,
       the plot will be saved in this file

    colors: list(str)
       A list of colors of mathplotlib

    """
            
    if len(experiment_list) > len(colors):
        error_msg =  "There can be at most {} plots.".format(len(colors))
        error_msg += "{} provided.".format(len(experiment_list))
        raise ValueError(error_msg)

    metagene_data = dump_metagene( ribo_handle     = ribo_handle, 
                                   site_type       = site_type,                           
                                   range_lower     = range_lower, 
                                   range_upper     = range_upper,
                                   sum_lengths     = True, 
                                   sum_references  = True,
                                   experiment_list = experiment_list )    
    
    normalization_factor = 1000000
    
    if normalize:
        ylabel = "Frequency (Per Million)"
        metagene_vectors = [ (metagene_data[lib].iloc[0] / get_total_reads(ribo_handle, lib) ) \
                                  * normalization_factor \
                             for lib in experiment_list]
    else:
        ylabel = "Frequency Raw"
        metagene_vectors = [ metagene_data[lib].iloc[0] for lib in experiment_list ]
    
    metagene_radius = get_metagene_radius(ribo_handle)
    plt.xticks(np.arange(-1*metagene_radius, metagene_radius+1, step=10 ) )

    for i, mg_vector in enumerate(metagene_vectors):
        obj1 = plt.plot(mg_vector,
               linestyle='-', color = colors[i] , linewidth = 1)
    
    plt.title(title)
    plt.xlabel("Relative Position")
    plt.ylabel(ylabel)
    if len(experiment_list) > 1:
        plt.legend(experiment_list, loc = (1.04,0))

    return metagene_data
    


@make_cli_function
def plot_metagene_wrapper( 
                   ribo_file, 
                   site_type, 
                   range_lower, range_upper,
                   experiment_list, 
                   output_file, 
                   title,
                   normalize    = False,
                   dump_to_file = "" ):
    """
    This is a wrapper function for plot_metagene
    See its documentation for details.

    This function is typically called by CLI
    and calls the actual plot function plot_metagene

    If dump_to_file is a non-empty string,
    the data used for the metagene plot is written
    in the provided file.
    """

    if not experiment_list:
        exit("Please provide experiment names.")


    with h5py.File(ribo_file, "r") as ribo_handle:
        if site_type.lower() == "start":
            site_type_actual = REF_DG_START_SITE_COV
        elif site_type.lower() == "stop":
            site_type_actual = REF_DG_STOP_SITE_COV
        else:
            raise ValueError("Wrong site type:", site_type)
        
        min_length, max_length = get_read_length_range(ribo_handle)

        plot_metagene( ribo_handle      = ribo_handle, 
                       site_type        = site_type_actual, 
                       range_lower      = range_lower, 
                       range_upper      = range_upper,
                       experiment_list = experiment_list, 
                       title            = title,
                       normalize        = normalize,
                       output_file      = output_file )
        

    if dump_to_file:
        dump_metagene_wrapper( ribo_file       = ribo_file, 
                               output_file     = dump_to_file, 
                               site_type       = site_type, 
                               sum_lengths     = True, 
                               sum_references  = True,
                               range_lower     = range_lower, 
                               range_upper     = range_upper,
                               experiment_list = experiment_list )

            
@fix_experiment_list
@save_plot
@check_experiment_list_in_ribo_handle
def plot_lengthdist( ribo_handle,
                     *, 
                     region_type, 
                     experiment_list, 
                     title, 
                     normalize   = False,
                     output_file = "",
                     colors      = PLOT_COLORS ):
    
    """
    Plots the distribution of the lengths of ribosome 
    footprints on a given region.

    If normalize is set to True, the values are divided
    by the sum and multiplied by 100 to give 
    the percentage on the y-axis.

    The number of experiments can not exceed the number
    of the elements os colors argument.

    Output file can be a pdf or a png file,

    Parameters
    ----------

    ribo_handle: h5py.File
      Open h5py handle

    region_type: str
      Transcript region to be plotted.
      It can be one of the UTR5, CDS, UTR3

    experiment_list: list(str)
      Name of the experiments to be plotted.

    title: str
      Title of the plot to appear at the top

    normalize: Boolean
       Determines whether raw frequencies 
       or percentages are to be plotted.

    output_file: str
       If provided, the resulting plot will be saved
       in this file. 
       pdf and png extensions are supported.

    colors: list(str)
       Colors of the individual experiments 
       to be passed on to matplotlib.

    """

    min_length, max_length = get_read_length_range(ribo_handle)

    region_data = dump_region(
                    ribo_handle     = ribo_handle, 
                    region_name     = region_type,
                    sum_lengths     = False, 
                    sum_references  = True, 
                    range_lower     = min_length , 
                    range_upper     = max_length ,
                    experiment_list = experiment_list)
    
    ylabel = "Frequency"
    
    if normalize:
        ylabel = "Frequency %"
        for name in experiment_list:
            length_dist       = np.array(region_data[name], dtype = np.int64)
            length_sum        = np.sum(length_dist)
            # Avoid division by zero
            length_sum        = 1 if length_sum == 0 else length_sum
            region_data[name] = ( length_dist / length_sum ) * 100

    ax = plt.subplot(111)
    plt.xticks(np.arange(min_length, max_length+1, step=2 ) )

    for i, name in enumerate(experiment_list):
        obj1 = plt.plot(region_data[name],
                   linestyle='-', color = colors[i] , linewidth = 1)
        
    if len(experiment_list) > 1:
        plt.legend(experiment_list, loc=(1.04,0))
    
    plt.title(title)
    plt.xlabel("Read Length")
    plt.ylabel(ylabel) 
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    return region_data


@make_cli_function
def plot_lengthdist_wrapper( 
                   ribo_file, 
                   region_type, 
                   experiment_list, 
                   output_file, 
                   title,
                   normalize,
                   dump_to_file = "" ):

    """
    This is a wrapper function for plot_lengthdist
    See its documentation for details.

    This function is typically called by CLI
    and calls the actual plot function plot_lengthdist

    If dump_to_file is a non-empty string,
    the data used for the metagene plot is written
    in the provided file.
    """
    
    with h5py.File(ribo_file, "r") as ribo_handle:
        min_length, max_length = get_read_length_range(ribo_handle)

        plot_lengthdist( ribo_handle      = ribo_handle, 
                         region_type      = region_type, 
                         experiment_list = experiment_list, 
                         title            = title, 
                         normalize        = normalize,
                         output_file      = output_file)


    if dump_to_file:
        dump_region_wrapper(ribo_file       = ribo_file, 
                            output_file     = dump_to_file, 
                            region_name     = region_type,
                            sum_lengths     = False, 
                            sum_references  = True, 
                            range_lower     = min_length, 
                            range_upper     = max_length,
                            experiment_list = experiment_list)


def _get_region_counts(   ribo_handle, 
                          range_lower, 
                          range_upper, 
                          experiment_list ):
    """
    Helper function for plotting region_counts
    
    It combines the data frames into one and fixes the keys.
    """
    region_getter = partial( dump_region, 
                             ribo_handle     = ribo_handle, 
                             sum_lengths     = True, 
                             sum_references  = True, 
                             range_lower     = range_lower, 
                             range_upper     = range_upper,
                             experiment_list = experiment_list)
                             
    region_list   = (UTR5_name, CDS_name, UTR3_name)

    region_counts = [ region_getter(region_name = r) for r in region_list ]
    
    #for i, region in enumerate(region_list):
    #    region_counts[i].index.name = region
        
    region_df     = pd.concat(region_counts, axis = 0, join = "outer")
    region_df.index = region_list
    
    return region_df    

@check_length_range
@fix_experiment_list
@save_plot
@check_experiment_list_in_ribo_handle
def plot_region_counts(ribo_handle, *,
                       experiment_list = [],
                       title            = "",
                       range_lower      = 0, 
                       range_upper      = 0, 
                       horizontal       = True,
                       output_file      = ""):
    """
    Generates barplots for the percentage of the reads
    of the regions UTR5, CDS and UTR3 for the provided
    experiments.

    Parameters
    ----------
    ribo_handle: h5py.File
      Open h5py handle

    experiment_list: list(str)
      Name of the experiments to be plotted.

    title: str
      Title of the plot to appear at the top

    output_file: str
       If provided, the resulting plot will be saved
       in this file. 
       pdf and png extensions are supported.

    Returns
    -------
    region_df : pd.DataFrame 
        Dataframe used to generate the plot.
        The actual counts (not the percentages).
    """

    region_df = _get_region_counts(   
                              ribo_handle     = ribo_handle, 
                              range_lower     = range_lower, 
                              range_upper     = range_upper, 
                              experiment_list = experiment_list )
    
    """
    region_getter = partial( dump_region, 
                             ribo_handle     = ribo_handle, 
                             sum_lengths     = True, 
                             sum_references  = True, 
                             range_lower     = range_lower, 
                             range_upper     = range_upper,
                             experiment_list = experiment_list)
    
    region_list   = (UTR5_name, CDS_name, UTR3_name)
    region_counts = [ region_getter(region_name = r) for r in region_list ]
    
    #for i, region in enumerate(region_list):
    #    region_counts[i].index.name = region
        
    region_df     = pd.concat(region_counts, axis = 0, join = "outer")
    region_df.index = region_list
    """
    percentage_df = region_df.divide(
                      other = region_df.sum(axis=0).T, axis = 1 ).T \
                                             * 100 
                                            
    cds_percentages = list(map("{0:.1f}".format, percentage_df[CDS_name]))    

    if horizontal:
        ax  = percentage_df.plot.barh(stacked=True)
            
        alignment = {'horizontalalignment': 'right', 
                     'verticalalignment'  : 'center'}
        for i, p in enumerate(cds_percentages):
            plt.text(50, i, p, **alignment)
    else:
        ax  = percentage_df.plot.bar(stacked=True)     
        
        alignment = {'horizontalalignment': 'center', 
                     'verticalalignment'  :   'baseline'}
        for i, p in enumerate(cds_percentages):
            plt.text(i, 50, p, **alignment)
        
    plt.legend(loc = (1.04,0))
    plt.title(title)
    
    return region_df.T
    


@make_cli_function
def plot_region_counts_wrapper(ribo_file, 
                               experiment_list, 
                               title,
                               range_lower, range_upper, 
                               output_file, 
                               dump_to_file,
                               horizontal = True):

    """
    A wrapper function for plot_region_counts

    This function is typically called by CLI.
    """

    with h5py.File(ribo_file, "r") as ribo_handle:
        count_df = plot_region_counts(
                           ribo_handle, 
                           experiment_list = experiment_list,
                           title            = title,
                           range_lower      = range_lower, 
                           range_upper      = range_upper, 
                           output_file      = output_file)

        if dump_to_file:
            region_df = _get_region_counts(   
                                      ribo_handle     = ribo_handle, 
                                      experiment_list = experiment_list,
                                      range_lower     = range_lower, 
                                      range_upper     = range_upper)     
            region_df.to_csv(dump_to_file, sep = DF_CSV_SEPARATOR)
