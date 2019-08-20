# -*- coding: utf-8 -*-
import functools
from functools import partial
from io import StringIO, IOBase

import matplotlib
import numpy as np
import pandas as pd
import h5py
import yaml

import matplotlib
import matplotlib.pyplot as plt
## Adjustment for larger plots in the jupyter notebook
plt.rcParams['figure.figsize'] = [15, 10]
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
# The "use" method is necessary for not seeing 
#gtk error message on the terminal
matplotlib.use('Agg')

from ..plot import plot_region_counts as main_plot_region_counts
from ..plot import plot_lengthdist as main_plot_lengthdist
from ..plot import plot_metagene as main_plot_metagene


from ..create import create_ribo
from ..settings import *
from ..info import ribo_info, get_experiment_info
from ..dump import dump_metagene, \
                   merge_metagene_dataframes, \
                   dump_region, \
                   dump_coverage

from ..core.exceptions import *
from ..core.get_gadgets import rnaseq_exists, \
                              has_coverage_data, \
                              get_experiment_names, \
                              has_coverage_data, \
                              rnaseq_exists, \
                              get_total_reads, \
                              get_reference_names, \
                              get_reference_lengths
from ..rnaseq import get_rnaseq
from ..metadata import get_metadata as main_get_metadata
from .api_helpers import _make_ribo_info_dict,\
                         _check_experiments
from .alias import ReferenceAlias

from .get_transcript_coverage import _get_transcript_coverage

                      

########################################################
###   A P I   F O R   R I B O G A D G E T S
########################################################

# Note in method definitions,
# always make the experiment (or experiments) parameters the SECOND parameter
# coming after self for consistency. Also this way we can easily
# decorate the functions with verifyer functions.  

class Ribo:
    """
    Ribo is an interface to ribo files.
    
    It provides access to ribo file attributes, metadata 
    and ribosome profiling data in a ribo file.
    
    Parameters
    ----------
    ribo_file: str, BytesIO
       Path to a ribo file or a handle to a ribo file
    file_mode: str , [choices: "r", "r+"]
       h5py file mode. Ribo file must exist. 
       So only read ("r") or read & write ("r+") modes are allowed.
       Be extremely careful with the r+ option.
       For most use cases, "r" (read-only) option should be sufficient. 
       You are strongly discouraged to delete 
       or modify existing attributes or data
       tables as this can corrupt the file.
    """
    
    def __init__(self, ribo_file, alias = None, file_mode = "r"):
        if file_mode not in ["r", "r+"]:
            raise RiboBaseError("Invalid file option."
                                "Options are r or r+.")
                                
        self._handle = h5py.File(ribo_file, file_mode)
            
        self._info             = _make_ribo_info_dict(self._handle)
        self._experiments      = get_experiment_names(self._handle)
        self._transcript_names = get_reference_names(self._handle)
        
        # Metadata is buffered in memory if user demands.
        self._metadata         = -1
        
        # Transcript indices, lengths and offsets are buffered in memory
        # Though this is not created on Ribo initialization
        # It is read from the file if access requested after object init.
        self._transcript_index   = {}
        self._transcript_lengths = {}
        self._transcipt_offsets  = {}
        
        if alias is None:
            self.alias = None
        else:
            self.alias = ReferenceAlias(
                              renaming_function = alias,
                              reference_names   = self._transcript_names)
        
    def __del__(self):
        self._handle.close()
        
    def close(self):
        self.handle.close()

    ######## P R O P E R T I E S ############################
        
    @property
    def experiments(self):
        """
        List of experiments in the ribo object
        """
        return self._experiments
    
    @property
    def reference_name(self):
        """
        Name of the reference (transacript assembly & annotation)
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_REFERENCE]]
    
    @property
    def format_version(self):
        """
        Ribo file format version
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_FORMAT_VERSION]]
    
    @property
    def ribopy_version(self):
        """
        Version of the ribopy used to creat the ribo file
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_VERSION]]
    
    @property
    def minimum_length(self):
        """
        Reads, used to generate the ribo file, are of at least minumum_length
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_LENGTH_MIN]]
    
    @property
    def maximum_length(self):
        """
        Reads, used to generate the ribo file, are of at least minumum_length
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_LENGTH_MAX]]
    
    @property
    def metagene_radius(self):
        """
        #nucleotides to the left and right of start / stop sites 
        in metagene coverage.    
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_METAGENE_RADIUS] ]
    
    @property
    def left_span(self):
        """
        #Nucleotides to the left of start / stop sites 
        to define UTR5 & UTR3 junction regions
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_LEFT_SPAN]]
    
    @property
    def right_span(self):
        """
        #Nucleotides to the right of start / stop sites 
        to define UTR5 & UTR3 junction regions
        """
        return self._info[RIBO_METADATA_FOR_DISPLAY[ATTRS_RIGHT_SPAN]]
        
    @property
    def info(self):
        """
        All Ribo attributes packed in a dictionary. 
        """
        return self._info
        
    @property
    def transcript_names(self):
        """
        Transcripts Names in the ribo reference 
        """
        return self._transcript_names
        
    @property
    def transcript_lengths(self):
        """
        Transcript lengths in dictionary form.
        
        transcript_name -> transcript_length
        """
        
        if not self._transcript_lengths:
            transcript_lengths       = get_reference_lengths(self._handle)
            self._transcript_lengths = dict( zip(self.transcript_names,
                                                 transcript_lengths) )
            
        return self._transcript_lengths
        
    @property
    def transcript_index(self):
        """
        Transcript indices in dictionary form.
        
        transcript_name -> transcript_length
        
        The transcript index coming from the order of Transcripts
        in the ribo file,
        """
        
        if not self._transcript_index:
            self._transcript_index = dict( map ( lambda x: (x[1], x[0]),
                                                 enumerate(self.transcript_names) ))
                                                    
        return self._transcript_index    
        
    @property
    def transcript_offsets(self):
        """
        Transcript offsets 
        
        This gives us the initial position of the transcript in determining coverage.
        If the transcripts are linearly lined up according to the order 
        in the ribo file, then this offset is the position of the first 
        nucleotide of the transcript.
        """
        offset = 0
        
        if not self._transcipt_offsets:
            for x in self.transcript_names:
                self._transcipt_offsets[x] = offset
                offset += self.transcript_lengths[x]
        
        return self._transcipt_offsets
        
    
    #/// P R O P E R T I E S   ////////////////////////
    
    def print_info(self, return_str = False):
        """
        Prints Ribo file information in string format.
        
        Parameters
        ----------
        return_str: bool
            If True, retuns the info_str and does not print.
            if False, prints the info_str and returns none.
        
        Returns
        -------
        info_str: str
            Ribo File Summary String
        
        Notes
        -----
        If this information is needed in a structured form,
        use the "info" attribute of the ribo object.
        """    
        info_str = ribo_info( self._handle )
        if return_str:
            return info_str
        else:
            print(info_str)
    
    #### D A T A   F R A M E    F U N C T I O N S ######## 
    
    def get_metagene(self,
                     site_type, 
                     experiments     = [], 
                     sum_lengths     = True, 
                     sum_references  = True, 
                     range_lower     = 0,
                     range_upper     = 0,
                     alias           = False 
                    ):
        """
        Returns metagene data at start / stop site
        
        Metagene data is reported for a range of read lengths.
        Every read length in this range is included in the returned data frame.
        This range is provided as an interval by specifying
        range_lower and range_upper, both of which are included.
        If range is not specified, the minimum and maximum 
        read lengths from the ribo file are taken as the range definition.
        
        Parameters
        ----------
        site_type: str [choices: 'start' , 'stop']
            Determines the site, around which metagene coverage is going to be
            reported. 
        sum_lengths: bool [default: True]
           If True, metagene data is summed up across the given
           length range.
           If False, metagene data is reported 
           for each length in the given range.
        sum_references: bool [default: False]
           If True, metagene data is summed across references (transcripts).
           If False, metagene data is reported for each transcript.
        experiments: list, str
           List of experiment(s). If empty, all experiments are included.
           
        Returns
        -------
        metagene_df: pandas.DataFrame
           This data frame contains coverage around start or stop site.
           The column labels are the relative positions, with respect
           to the start  or stop site.
           The index of the dataframe depends on the two parameters:
           sum_lengths and sum_references.
           Whichever of these parameters are set to True, won't be in the index.
           Therefore, if sum_lengths is True, there won't be a length index
           in the metagene_df because the values are summed across the lengths. 
        """
        
        if type(experiments) == str:
            experiment_list = [experiments]
        else:
            experiment_list = experiments
        
        internal_site_type = ""
        
        if site_type.lower() == "start":
            internal_site_type = "start_site_coverage"
        elif site_type.lower() == "stop":
            internal_site_type = "stop_site_coverage"
        else:
            raise InvalidName("Site type can be either \"start\" or \"stop\".")
        
        if alias:
            if self.alias is None:
                raise AliasError("This ribo object doesn't have a defined alias.")
            alias_list = self.alias.aliases
        else:
            alias_list = None
            
        metagene_df = \
               dump_metagene( ribo_handle     = self._handle, 
                              site_type       = internal_site_type, 
                              sum_lengths     = sum_lengths, 
                              sum_references  = sum_references,
                              range_lower     = range_lower, 
                              range_upper     = range_upper,
                              experiment_list = experiment_list,
                              alias           = alias_list )
        
        metagene_df = merge_metagene_dataframes(
                                  df_dict        = metagene_df, 
                                  sum_references = sum_references, 
                                  sum_lengths    = sum_lengths)
             
        return metagene_df   

    
    def get_region_counts(self, 
                          region_name, 
                          sum_lengths     = True, 
                          sum_references  = True, 
                          range_lower     = 0,
                          range_upper     = 0,
                          experiments     = [],
                          alias           = False):
        """
        Returns number of reads mapping to UTRs or CDS. 
        
        Region counts are reported for a range of read lengths.
        Every read length in this range is included in the returned data frame.
        This range is provided as an interval by specifying
        range_lower and range_upper, both of which are included.
        If range is not specified, the minimum and maximum 
        read lengths from the ribo file are taken as the range definition.
        
        Parameters
        ----------
        region_name: str [choices: "UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3"]
           For the definition of these regions, check the documentation at
           https://ribopy.readthedocs.io
           The region definitions are coming from the annotation in the ribo file.
           
        sum_lengths: int [default: True]
           If True, sum region counts across read lengths.
           If False, report region counts for each length.
        sum_references: int [default: True]
           If True, sum region counts across references (transcripts).
           If False, report region counts for each reference.
        range_lower: int
           Minimum read length to be included in the result.
         range_upper: int
           Maximum read length to be included in the result.
         experiments: list, str [default: []]
           List of experiment(s). If empty, all experiments are included.
           
        Returns
        -------
        region_counts: pd.DataFrame
           The index of the dataframe depends on the parameters
           sum_lengths and sum_references. More precisely, if sum_lengths
           is False, there will be an index for read length.
           If True, there won't be such an index as data is aggregated by
           read lengths. Similarly, if sum_references is False, there will be
           an index for reference (transcript) names. If True, there won't be
           such an index.
           The columns correspond to experiments.
           
        """
                          
        if type(experiments) == str:
            experiments = [experiments]
            
        if alias:
            if self.alias is None:
                raise AliasError("This ribo object doesn't have a defined alias.")
            alias_list = self.alias.aliases
        else:
            alias_list = None
        
        region_counts = dump_region(ribo_handle      = self._handle,
                                    region_name      = region_name,
                                    sum_lengths      = sum_lengths, 
                                    sum_references   = sum_references, 
                                    range_lower      = range_lower , 
                                    range_upper      = range_upper ,
                                    experiment_list  = experiments,
                                    alias            = alias_list)
        return region_counts
        
        
    def get_length_dist(
                    self,    
                    region_name,
                    experiments  = []):
        """
        Returns the number of reads for each length for a given region
        
        Parameters
        ----------
        region_name : str, list
           
        """
        if region_name not in EXTENDED_REGION_names:
            error_message = "Invalid region name: {}\n".format(region_name) + \
                            "The region name should be one of the following:\n" + \
                            ", ".join(EXTENDED_REGION_names) 
            raise InvalidName(error_message)            
                    
        length_dist = dump_region(
                        ribo_handle      = self._handle,
                        experiment_list  = [],
                        region_name      = region_name,
                        sum_lengths      = False, 
                        sum_references   = True, 
                        range_lower      = self.minimum_length , 
                        range_upper      = self.maximum_length )
                        
        return length_dist
                        
    
    #////// D A T A   F R A M E    F U N C T I O N S   ///////////////////
    
    ######   P L O T  F U N C T I O N S   ######################
    
    def plot_metagene( self, 
                       site_type, 
                       title           = "",
                       experiments     = [], 
                       range_lower     = 0, 
                       range_upper     = 0,
                       normalize       = False,
                       output_file     = "",
                       colors          = PLOT_COLORS ):
        """Generates coverage plots around start / stop sites.
        
        "metagene_radius" many nucleotides are taken on either side of the
        site type. The axis is the positions of the nucleotides relative to
        "site_type" and the y-axis is their frequency.
        
        Setting normalize = True might be helpful when plotting 
        more than one experiment. Values are normalized by total number 
        of mapped reads in the experiment, 
        
        Arguments
        ---------
        site_type: str [choices: "start", "stop"]
           Coverage plot is centered around the site type.
        title: str, default = ""
           Title of the plot       
        experiments: list, str
            List of experiments to be plotted        
        range_lower: int
            Minimum read length to be included in the metagene coverage data.        
        range_upper: int
            Maximum read length to be included in the metagene coverage data.        
        normalize: bool (default = False)
            Normalize metagene data by number of mapped reads in the experiment     
        output_file: str, (default = "")
            If non-empty, the plot will be saved in the provided path.    
        colors: list, (default: ["blue",  "red",    "green", "brown", "orange", "violet", "black"])
            Colors of lines.
        """
        
        if type(experiments) == str:
            experiments = [experiments]
            
        internal_site_type = ""
        
        if site_type.lower() == "start":
            internal_site_type = "start_site_coverage"
        elif site_type.lower() == "stop":
            internal_site_type = "stop_site_coverage"
        else:
            raise InvalidName("Site type can be either \"start\" or \"stop\".")
        
        return  main_plot_metagene( 
                                        self._handle, 
                                        site_type       = internal_site_type, 
                                        experiment_list = experiments, 
                                        title           = title,
                                        range_lower     = range_lower, 
                                        range_upper     = range_upper,
                                        normalize       = normalize,
                                        output_file     = output_file,
                                        colors          = PLOT_COLORS )
                                        

    def plot_lengthdist( self, 
                         region_type, 
                         experiments, 
                         title       = "", 
                         normalize   = False,
                         output_file = "",
                         colors      = PLOT_COLORS ):
        """Generates distribution of the reads according to length
        
        The x-axis is the read length and the y-axis is the number of reads
        mapping to the particular region.
        
        Parameters
        ----------
        region_name: str [choices: "UTR5", "UTR5_junction", "CDS", "UTR3_junction", "UTR3"]
            For the definition of these regions, check the documentation at
            https://ribopy.readthedocs.io
            The region definitions are coming from the annotation in the ribo file.
        experiments: list, str
            List of experiments to be plotted
        title: str
            Title of the plot
        normalize: bool (default = False)
            Normalize each experiment by total number of reads.
        output_file: str (default = "")
            If provided, the output will be saved in this path.
        colors: (default: ["blue",  "red",    "green", "brown", "orange", "violet", "black"])
            Colors of lines.
        """
        
        return  main_plot_lengthdist( 
                                ribo_handle     = self._handle,
                                region_type     = region_type, 
                                experiment_list = experiments, 
                                title           = title, 
                                normalize       = normalize,
                                output_file     = "",
                                colors          = PLOT_COLORS )
    
    def plot_region_counts(
                       self,
                       experiments,
                       title            = "",
                       range_lower      = 0, 
                       range_upper      = 0, 
                       horizontal       = True,
                       output_file      = ""):
        """Generates bar plots of region counts
        
        The bar plot coming from the percentages of the counts of the regions:
        UTR5, CDS and UTR3
        
        Parameters
        ----------
        experiments: list, str
            List of experiments to be plotted
        title: str
            Title of the plot
        range_lower: int
            Minimum read length to be included in the bar plot data.        
        range_upper: int
            Maximum read length to be included in the bar plot data.  
        horizontal: bool (default = True)
            Generates bar plots horizontally. 
            Especially for long experiment names,  this is the preferred method.
        output_file: str (default = "")
            If provided, the output will be saved in this path.
        """
        
        return main_plot_region_counts(
                                  ribo_handle      = self._handle,
                                  experiment_list  = experiments,
                                  title            = title,
                                  range_lower      = range_lower, 
                                  range_upper      = range_upper, 
                                  horizontal       = horizontal,
                                  output_file      = output_file)

    #/////   P L O T  F U N C T I O N S   //////////////////////
    
    #### R N A S E Q   F U N C T I O N S   ######################
    
    def has_rnaseq(self, experiment):
        """
        Does experiment has RNA-Seq data?
        
        Parameters
        ----------
        experiment: str
           Name of the experiment
        Returns
        -------
        rnaseq_exists: bool
        """
        return self._info["experiments"][experiment]["RNA-Seq"]
    
    def get_rnaseq(self, experiments = None):
        """
        Returns region counts coming from RNA-Seq data.
        
        Note that RNA-Seq data is an optional entity of Ribo File. 
        
        In contrast to region counts for ribosome profiling data,
        the resulting data frame do not have separate entries for read lengths.
        The columns of the data frame correspond to regions.
        The index of the Data Frame has two levels: experiment, transcript. 
        
        parameters
        ----------
        experiments: list, str
           List of experiment(s) whose RNA-Seq data is to be reported
        rnaseq_df: pd.DataFrame
           Data Frame containing RNA-Seq counts for each region.
        """
        
        rnaseq_df =  get_rnaseq(self._handle, experiments)
        return rnaseq_df
    
    ######## M E T A D A T A   F U C N T I O N S  ######################
    @_check_experiments
    def has_metadata(self, experiment = None):
        """
        Checks if Ribo has metadata
        
        Parameters
        ----------
        experiment: str (default = None)
           If None, the metadata of the ribo file itself is inquired.
           Otherwise, checks if the given experiment has metadata
        Returns
        -------
        result: bool
        """
        if experiment:
            result = bool(self.info[EXPERIMENTS_name][experiment]["Metadata"])
        else:
            if self._metadata != -1 and not self._metadata:
                result = False
            else:
                result = True
            
        return result
        
    def get_metadata(self, experiment = None):
        """
        Returns user defined metadata in dictionary form
        
        Parameters
        ----------
        experiment: str (default = None)
           If None, the metadata of the ribo file itself is returned.
           Otherwise, returns the metadata of the given experiment.
        Returns
        -------
        result: dict
        """
        
        if experiment:
            metadata = main_get_metadata(self._handle, name = experiment)
        else:
            # A minor optimization
            # When the ribo file is created, the metadata is NOT read from disk
            # But upon the first inquiry to the metadata of the root ribo
            # file, we save it in an object attribute called _metadata.
            # Thus root ribo metadata is read from disk only once.
            # For experiments, metadata is read from disk every time. 
            if self._metadata == -1:
                self._metadata = yaml.load(self._handle.attrs.get(USER_METADATA, ""), 
                                           Loader=yaml.FullLoader)
                                     
            metadata = self._metadata
        return metadata
                
        
    #/////// M E T A D A T A   F U C N T I O N S  ######################
    
    
    ### C O V E R A G E    F U N C T I O N S ##########################
    @_check_experiments
    def has_coverage(self, experiment):
        """
        return has_coverage_data(self._handle, experiment_name)
        
        Parameters
        ----------
        experiment: str
            Named of the experiment whose coverage data is inquired.
        Returns
        -------
        result: bool
        """

        return self.info[EXPERIMENTS_name][experiment]["Coverage"]
    
    def get_coverage( self, 
                      experiment,
                      range_lower = 0, 
                      range_upper = 0,
                      alias       = False):
        """
        Returns coverage at nucleotide resolution.
        
        Note that RNA-Seq data is an optional entity of Ribo File.
        
        """
        
        if alias:
            if self.alias is None:
                raise AliasError("This ribo object doesn't have a defined alias.")
            alias_list = self.alias.aliases
        else:
            alias_list = None
        
        return dump_coverage( ribo_handle     = self._handle,
                              experiment_name = experiment,
                              range_lower     = range_lower, 
                              range_upper     = range_upper,
                              alias           = alias_list)
    
    @_check_experiments                          
    def get_transcript_coverage( self, 
                                 experiment,
                                 transcript, 
                                 range_upper = 0,
                                 range_lower = 0,
                                 sum_lengths = False,
                                 alias       = False):
                                 
        if alias:
            if self.alias is None:
                raise AliasError("This ribo object doesn't have a defined alias.")
            transcript_name = self.alias.get_original_name(transcript)
        else:
            transcript_name = transcript
                                                
        #transcript_index        = self.transcript_index.get(transcript, None)
        transcript_length       = self.transcript_lengths.get(transcript_name, None)
        
        if transcript_length == None:
            raise(TranscriptDoesntExist( "{} not found".format(transcript) ))
        total_transcript_length = np.sum( tuple(self.transcript_lengths.values()))                 
                                 
        transcript_coverage = _get_transcript_coverage(
                                     ribo              = self,
                                     experiment_name   = experiment,
                                     transcript        = transcript_name,
                                     transcript_length = transcript_length,
                                     total_length      = total_transcript_length,
                                     range_lower       = range_lower, 
                                     range_upper       = range_upper,
                                     sum_lengths       = sum_lengths)
                                    
        return transcript_coverage
