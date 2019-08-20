# -*- coding: utf-8 -*-
from collections import OrderedDict
import time
import os
import re

import h5py
import pandas as pd
import numpy as np

from .settings import *
from .io.file import get_alignment_file_handle
from .core.get_gadgets import *
from .core.metagene import *
from .metadata import set_metadata
from .core.exceptions import *
from .core.create_experiment import create_experiment
from .core.verify import make_cli_function
from .io.file import (open_by_extension,
                      read_all_lines, 
                      flex_open) 
from ._version import __format_version__,\
                      __version__ 

########################################################


def _get_experiment_name_from_filepath(filepath):
    """
    Extracts name from the ribo file path
    after stripping the file extension
    """

    basename        = os.path.basename(filepath)
    basename_pieces = basename.split(".")
    experiment_name = ".".join(basename_pieces[:-1])
    
    return experiment_name

def _check_experiment_name(name):
    """
    Raises an exception if the name is invalid.
    """
    
    error_message = \
    """
    Invalid Experiment Name! 
    Name must contain at least one alphanumeric character and
    it can only contain alphanumeric characters, '_', '-' and '.'.
    """
    
    alphanum_test   = re.match("[a-zA-Z0-9]+", name)
    valid_char_test = re.fullmatch("[a-zA-Z0-9\_\-.]+", name)
    
    if alphanum_test and valid_char_test:
        return True
    else:
        raise InvalidName(error_message)


def initialize(h5_handle, experiment_name, reference_name):
    """
    Basic initialization of a ribo file.
    It creates the hdf5 groups for experiments
    and reference.  

    Parameters
    ----------
    h5_handle: h5py.File
        Handle for writable hdf5 file

    experiment_name : str 
        A file group under
        experiments is created with this name

    reference_name : str 
        User defined name of the reference

    Returns
    -------
    None

    """

    _check_experiment_name(experiment_name)
    
    h5_handle.create_group(EXPERIMENTS_name)
    h5_handle.create_group(REFERENCE_name)
    h5_handle[REFERENCE_name].attrs[ATTRS_REFERENCE]      = reference_name
    h5_handle[REFERENCE_name].attrs[ATTRS_FORMAT_VERSION] = __format_version__
    h5_handle[REFERENCE_name].attrs[ATTRS_VERSION]        = __version__

    exp_group = h5_handle[EXPERIMENTS_name].create_group(experiment_name)



def set_reference_names_and_lengths(h5_handle, lengths_file):
    """
    Creates reference names and lengths fields in the references part
    of the ribo file.

    The order (of transcripts) is crucial because 
    the same order will be kept in all parts of the file.  


    Parameters
    ----------

    h5_handle: h5py.File
        Handle for writable hdf5 file

    lengths_file : str 
        Tab separated file of two columns
        First column: Transcript names
        Second Column: Transcript lengths
        One line per transcript.

    Returns
    -------
    Reference names and lengths
    """

    reference_handle = h5_handle[REFERENCE_name]
    reference_df     = pd.read_csv(lengths_file, sep = "\s+", 
                               names = ["name","length"])
    name_lengths     = list(map(len, reference_df["name"]) )
    max_len          = max(name_lengths)
    dt               = np.dtype("S" + str(max_len))
    ref_names        = np.array( reference_df["name"], dtype = dt )
    ref_lengths      = np.array( reference_df["length"], 
                                 dtype = np.uint32 )

    reference_handle.create_dataset( 
        REF_DG_REFERENCE_NAMES, 
        data        = ref_names, 
        dtype       = dt, 
        compression = DEFAULT_COMPRESSION, 
        fletcher32  = DEFAULT_FLETCHER32 )

    reference_handle.create_dataset( 
        REF_DG_REFERENCE_LENGTHS, 
        data        = ref_lengths, 
        dtype       = np.uint32, 
        compression = DEFAULT_COMPRESSION, 
        fletcher32  = DEFAULT_FLETCHER32 )
    
    return ( get_reference_names(h5_handle) , 
             get_reference_lengths(h5_handle))

#########################################################

def post_check_annotation(h5_handle, pre_annotation):
    """
    More sanity checks for the annotation.
    See the function "check_annotation" for details.
    """ 

    error_base     = "Annotation error:\n"
    error_message  = ""
    ref_lengths    = get_reference_lengths(h5_handle)
    ref_names      = get_reference_names(h5_handle)
    len_dict       = dict( tuple(zip( ref_names, ref_lengths )) )

    for this_name, boundaries in pre_annotation.items():
        utr5, cds , utr3 = tuple(boundaries.values())

        if cds[0] == cds[1]:
            error_message +=  this_name + \
                              " : must have a positive length CDS\n"

        if utr3[0] == utr3[1]:
            error_message +=  this_name + \
                              " : must have a positive length UTR3\n"

        if utr5[1] != cds[0]:
            error_message += this_name + \
                " : has gaps or overlaps in annotation. Check its boundaries.\n"

        if cds[1] != utr3[0]:
            error_message += this_name + \
                " : has gaps or overlaps in annotation. Check its boundaries.\n"

        if utr3[1] != len_dict[this_name]:
            error_message += this_name + \
                " : utr3 3 prime (right) boundary " + \
                "must be equal to the gene length.\n"

        if error_message:
            raise AnnotationError(error_base + error_message)





def check_annotation( h5_handle, annotation_lines ):
    """
    Performs some sanity checks in the given annotation file.

    Note that the transcript (reference) names and lengths
    are provided in a tsv file and recorded prior to this
    step.

    Every trasnscript must have an annotation
    and a non-zero length of CDS and UTR3 regions.

    Also there shouldn't be any gaps in the annotation.
    In other words, there shouldn't be a region of a
    transcript that bnelongs to none of CDS, UTR5 or UTR3.

    End position of one region must be equal to the start
    position of the next region. Note that the region
    boundaries are in bed form. This means, it is 0-based
    and start position is inclusive and the end position is exclusive.

    It raises an AnnotationError once it detects an inconsistency
    for a transcript.

    Parameters
    ----------

    h5_handle:  h5py.File 
        Handle to an open, readable/ writable hdf5 file

    annotation_lines: list
        An array of bed file entries

    """
    ref_lengths = get_reference_lengths(h5_handle)
    ref_names   = get_reference_names(h5_handle)

    line_number   = 0
    error_message = "Annotation error in line number: "

    recorded_annotations = OrderedDict()
    for r in ref_names:
        recorded_annotations[r] = \
            OrderedDict( ( (UTR5_name, [0,0]), 
                         (CDS_name, [0, 0]), 
                         (UTR3_name, [0,0]) ) )


    for this_line in annotation_lines:
        line_number +=1
        contents     = this_line.strip().split()
        if len(contents) < 6:
            continue

        current_ref = contents[0]
        start       = int(contents[1])
        end         = int(contents[2])
        region      = contents[3].upper()
        strand      = contents[5]

        if region not in REGION_names:
            error_message += str(line_number) + "\n"
            error_message += "\nInvalid region identifier:\n"
            error_message += "{} not in {}".format(region, REGION_names)
            raise AnnotationError(error_message)

        if start > end:
            error_message += str(line_number) + "\n"
            error_message += "\nInvalid region boundary:\n"
            error_message += "{} > {}".format(start, end)
            raise AnnotationError(error_message)

        if current_ref not in ref_names:
            error_message += str(line_number) + "\n"
            error_message += "Invalid reference name {}".\
                                format(current_ref)
            raise AnnotationError(error_message)


        recorded_annotations[current_ref][region] = [start, end]

    post_check_annotation( h5_handle, recorded_annotations)



def set_annotation( h5_handle, annotation_lines ):
    """
    Note that everything is in the transcriptomic coordinates
    and everything is with respect to + strand

    In the annotation we store a matrix of dimensions N * 3
    Each row is for one reference (in the same order as
    reference names and lengths)
    The entries are coming from BED coordinates
    <UTR5 NED POS> , <CDS END POS>, <UTR3 END POSITION>
    Since UTR5 begins with 0, and the end of one region is
    the start of the next one, the above numbers give us all
    region boundaries. 

    Parameters
    ----------

    h5_handle: h5py.File
        Handle to an open, readable/ writable hdf5 file

    annotation_lines: list 
        An array of bed file entries

    Returns
    -------

    boundaries: list( tuple )

                An array of triplets. Each triplet corresponds
                to the UTR5, CDS and UTR3 boundaries of a transcript.
                The order of the transcripts are the same as
                their order in get_reference_names(ribo)

                Thus each triplet is of the form
                ( (UTR5_left, UTR5_right), (CDS_left, CDS_right), 
                  (UTR3_left, UTR3_right) )

                These coordinates are in bed conventions.
                Therefore the left coordinates ARE inclusive 
                and the right coordinates are NOT.
                The coordinates are 0-based.

                Also,

                UTR5_right = CDS_left
                CDS_right  = UTR3_left
                UTR3_right = ref_length

                must hold.
    """

    # Raises an AnnotationError exception if 
    # there is a porblem with the annotation 
    check_annotation(h5_handle, annotation_lines)

    number_of_refs = get_number_of_references(h5_handle)
    region_dict    = OrderedDict()

    # Initialize all regions to zeros
    for region in REGION_names:
        region_dict[region] = np.full(shape      = (number_of_refs), 
                                      fill_value = 0, 
                                      dtype      = np.uint32)

    ref_index   = 0
    ref_list    = get_reference_names(h5_handle)
    current_ref = ref_list[0]
    prev_ref    = current_ref


    for this_line in annotation_lines:
        contents = this_line.strip().split()
        if len(contents) < 6:
            continue
        current_ref = contents[0]
        start       = int(contents[1])
        end         = int(contents[2])
        region      = contents[3].upper()
        strand      = contents[5]

        # increment the index if we see a new reference
        # make sure that the next ref in the bed file is the same as the\
        # ref in the ref_list
        # this also ensures that the bed file is in the correct order.

        if current_ref != prev_ref:
            ref_index +=1
            assert current_ref == ref_list[ref_index]
            prev_ref = current_ref

        region_dict[region][ref_index] = (end)


    concatenated_arrays = np.concatenate( tuple(region_dict.values()) )
    ref_array           = np.reshape(concatenated_arrays, 
                                      newshape = (number_of_refs, 3), 
                                      order    = "F" )

    h5_handle[REFERENCE_name].create_dataset(
                                REF_ANNOTATION_NAME,
                                data        = ref_array,
                                dtype       = np.uint32,
                                compression = DEFAULT_COMPRESSION, 
                                fletcher32  = DEFAULT_FLETCHER32)

    boundaries = tuple(map( lambda x: 
                              [(0, x[0]) , (x[0], x[1]), (x[1], x[2])] , 
                            ref_array ) )

    return boundaries
  
##########################################################

###########################################################    

def set_coverage_vectors_individual(pivots, radius):
    """
    Computes the total contribution of all transcripts
    around start sites.
    The pivots iterate over transcripts.
    Note that it may well happen that  some transcripts
    are short on the 5' or 3' ends so that
    they may not contribute to all radius many nucleotides on
    any end of the start / stop site.
    So we initialize contribution to 0, and increement the
    parts that are covered by the transcript.
    
    Example:
    If our transcript T is as follows:
    (Bed format coordinates)
    UTR5: 0 2
    CDS: 2 10
    UTR3: 10 15

    if radius is 3,
    its corresponding contribution vector will be
    w = (0,1,1,1,1,1,1)
    Note that the left-most contribution is 0
    because this transcript has only 2 nucleotides on UTR5

    This vector 2 will be added to result, 
    along with contribution of other transcripts,
    to compute the overal contribution around the start site.  

    """
    result = np.zeros( 2*radius + 1, dtype = np.uint32 )

    for p in tuple(pivots):
        this_coverage = np.zeros( 2*radius + 1, dtype = np.uint32 )
        start_index   = max( radius - p[1] , 0)
        end_index     = radius + 1 + min( radius , p[2] - p[1] )

        this_coverage[start_index:end_index] = 1
        result += this_coverage

    return result

###########################################################    

def set_coverage_vectors(ribo, radius):
    """
    Coverage vectors hold the contribution of the transcripts
    around start / stop sites .
    These vectors can be used to normalize metagene data
    ( start / stop sites) later. 
    """

    boundaries        = get_region_boundaries(ribo)
    start_site_pivots = map( lambda x: ( 0, x[1][0], x[2][1] ) , 
                               boundaries )
    stop_site_pivots  = map( lambda x: ( 0, x[1][1], x[2][1] ) , 
                               boundaries )

    start_site_boundaries = \
       set_coverage_vectors_individual(start_site_pivots, radius)
    stop_site_boundaries = \
       set_coverage_vectors_individual(stop_site_pivots, radius)

    ribo[REFERENCE_name].create_dataset(
                           REF_DG_START_SITE_COV, 
                           data       = start_site_boundaries, 
                           dtype      = np.uint32,
                           fletcher32 = DEFAULT_FLETCHER32 )
    ribo[REFERENCE_name].create_dataset(
                           REF_DG_STOP_SITE_COV, 
                           data       = stop_site_boundaries, 
                           dtype      = np.uint32,
                           fletcher32 = DEFAULT_FLETCHER32 )

###########################################################

def set_ribo_info(ribo, 
                  length_min, length_max, 
                  reference_name,
                  metagene_radius, 
                  left_span, right_span):

    """
    Essential ribo metadata is stored 
    in the atrributes of the hdf5 file.
    These attributes are crucial and shouldn't be changed
    after file creation.
    """

    ribo.attrs[ATTRS_VERSION]         = __version__
    ribo.attrs[ATTRS_FORMAT_VERSION]  = __format_version__

    ribo.attrs[ATTRS_LENGTH_MIN]      = length_min
    ribo.attrs[ATTRS_LENGTH_MAX]      = length_max
    ribo.attrs[ATTRS_TIME]            = time.time()
    ribo.attrs[ATTRS_METAGENE_RADIUS] = metagene_radius
    ribo.attrs[ATTRS_LEFT_SPAN]       = left_span
    ribo.attrs[ATTRS_RIGHT_SPAN]      = right_span
    ribo.attrs[ATTRS_REFERENCE]       = reference_name

##########################################################

def _get_metadata(file_path):
    if not file_path:
        return ""
    else:
        with flex_open(file_path) as metastream:
            meta_str = metastream.read()
        return meta_str

##########################################################    

def create_ribo(ribo, 
                experiment_name, 
                alignment_file,
                reference_name,
                lengths_file, 
                annotation_file,
                metagene_radius, 
                left_span, right_span,
                length_min, length_max,
                alignment_format    = 'bed',
                store_coverage      = False,
                ribo_metadata       = None,
                experiment_metadata = None,
                nprocess            = 1,
                tmp_file_prefix     = ""):
    """
    Creates a ribo file for a single experiment.


    The essential input data of this function are
       ribo file handle, annotation, transcript length file,
       alignment file

    The essential input, as single value parameters, are:
       experiment name and reference name,
       metagene_radius, left_span, right_span,
       length_min, length_max.    

    It writes annotation, metadata and 
    ribosome profiling data of the experiment
    to the given ribo file handle.

    Parameters
    ----------
    ribo :  h5py.File 
       Handle for the open hdf5 file

    experiment_name : str
               The name of the single experiment in the ribo file.
               An datagroup under experiments will be created accordingly 
    
    lengths_file : str
               File path. 
               A tab separated file containing the reference lengths
               The first column holds the  reference names
               The second column holds the reference lengths
               Note that the reference names 
               should be exactly the same as in
               the alignment file and the annotation file
               in the SAME ORDER

    alignment_file : str
                 File path.
                 Can be a bed, bam or sam file.
                 It can also be an StringIO stream.
                 If this is an empty string,
                 the input will be read from the standard io


    reference_name : str
             User defined name for the reference set:
             (transcript sequence, annotation)

    annotation_file : str
                  File path.
                  A bed file annotation the UTR5, CDS and UTR3
                  regions of the transcirpts. 
                  This can also be a StringIO object

    metagene_radius : int
                  The number of nucleotides to the left and right of
                  the start / stop sites in the metagene analysis.

    left_span : int
            The number of nucleotides to the left of start / stop sites
            for defining UTR5 & UTR3 junction regions.

    right_span : int
             The number of nucleotides to the right of start / stop sites
             for defining UTR5 & UTR3 junction regions.

    length_min : int
             Minumum read length to be counted in the ribosome profiling.
             Inclusive.

    length_max : int
             Maximum read length to be counted in the ribosome profiling.
             Inclusive. Note that all read lengths from 
             length_min to lenght max are counted.

    store_coverage: Boolean
        Store the coverage, at length level, for all transcripts. 

    ribo_metadata : str, dict or IOStream
        User-provided metadata of the ribo file in yaml format.

    experiment_metadata : str, dict or IOStream
        User-provided metadata of the experiment in yaml format.

    nprocess : int
               Number of parallel processes for 
               extracting the ribosome profiling data.

    tmp_file_prefix : int
                      File name prefix for the intermediate files
                      created (and later deleted).
                      For example, reads of the same length 
                      are collected in one inetermediate file. 

    Returns
    -------
    None

    """

    initialize(ribo, experiment_name, reference_name)

    (ref_names, ref_lengths) = \
       set_reference_names_and_lengths(ribo , lengths_file)

    annotation = set_annotation( ribo, 
                                 read_all_lines(annotation_file) )

    set_coverage_vectors(ribo, metagene_radius)

    set_ribo_info(ribo, 
                  length_min, length_max,
                  reference_name, 
                  metagene_radius, 
                  left_span, right_span)

    alignment_file_handle = \
            get_alignment_file_handle(alignment_file, 
                                      alignment_format = alignment_format )

    set_metadata(ribo_handle = ribo, 
                 name        = "", 
                 metadata    = ribo_metadata)

    create_experiment(
        ribo_exp_handle       = ribo[EXPERIMENTS_name][experiment_name], 
        experiment_name       = experiment_name, 
        alignment_file_handle = alignment_file_handle,
        metagene_radius       = metagene_radius,
        ref_names             = ref_names,
        ref_lengths           = ref_lengths,
        region_coordinates    = annotation ,
        left_span             = left_span, 
        right_span            = right_span,
        length_min            = length_min, 
        length_max            = length_max,
        store_coverage        = store_coverage,
        metadata              = experiment_metadata, 
        nprocess              = nprocess,
        tmp_file_prefix       = tmp_file_prefix)


        

############################################################################

def create_ribo_file(ribo_file_path, experiment_name, 
                     *argv, **kwargs):
    """
    Create a ribo file in the given file path.
    This is a wrapper for the create_ribo function to be used by cli.
    See create_ribo function for details.

    Parameters
    ----------

    ribo_file_path: str or path
        The path of the ribo file

    experiment_name: str
        The name of the experiment in the ribo file

    *argv, **kwargs: passed on to create_ribo function

    """

    kwargs["ribo_metadata"]       = \
         _get_metadata(kwargs.get("ribo_metadata"))
    kwargs["experiment_metadata"] = \
         _get_metadata(kwargs.get("experiment_metadata"))
    tmp_file_prefix = ribo_file_path + "."

    if not experiment_name:
        experiment_name = _get_experiment_name_from_filepath(ribo_file_path)


    with h5py.File( ribo_file_path, "w" ) as ribo_handle:
        try:
            create_ribo(ribo_handle,
                        experiment_name,  
                        *argv, **kwargs, 
                        tmp_file_prefix = tmp_file_prefix)
        except RiboBaseError as e:
            ribo_handle.close()
            os.remove(ribo_file_path)
            raise RiboBaseError(e)
        except:
            print("An unexpected error happend!")
            ribo_handle.close()
            os.remove(ribo_file_path)
            raise
            

############################################################################

@make_cli_function
def create_ribo_file_wrapper(*argv, **kwargs):
    """
    Wrapper function for create_ribo_file
    to be used by CLI
    """

    try:
        create_ribo_file(*argv, **kwargs)
    except RiboBaseError as e:
        print("Ribo file creation has failed!")
        print(e)
        return 1

    return 0
