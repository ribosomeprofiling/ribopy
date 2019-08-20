# -*- coding: utf-8 -*-

import gzip
import io
from io import StringIO, IOBase
import sys
import pipes
import shutil

from ..core.exceptions import RiboBaseError

#############################################

def get_file_opener_by_extension(file_path):
    file_string = file_path.lower()
    if file_string.endswith("gz") or \
       file_string.endswith("gzip"):

       return gzip.open
    else:
        return open

#####################################################################        

def open_by_extension(file_path, mode = "rt"):
    file_opener = get_file_opener_by_extension(file_path)
    return file_opener(file_path, mode)

#####################################################################

def read_all_lines(file_input):
    if issubclass( type(file_input), IOBase ):
        return file_input.readlines()

    elif type(file_input) == str:
        with open_by_extension( file_input ) as file_handle:
            all_lines = file_handle.readlines()
        return all_lines
    
    else:
        raise RiboBaseError("Input file must be a stream or a file path")

#####################################################################

def get_alignment_file_handle(alignment_file, alignment_format = "bed"):

    # Probably we dont need this
    # but an extra check doesnt hurt
    if alignment_format not in ("bed", "bam"):
        raise RiboBaseError("File format can be bam or bed.")

    # If the alignment file is bam
    # it needs to be converted to bed using bamToBed
    if alignment_format == "bam":
        if not shutil.which('bamToBed'):
            raise RiboBaseError("RiboPy could not "
                                "find the executable bamToBed."
                                "This executable is needed for bam files.")
        pipe_t = pipes.Template()

    if not alignment_file:
        # ALignment data is coming from standard input
        if alignment_format == "bed":
            alignment_file_handle = sys.stdin
        else:
            # We need to convert bam to bed using
            # bamToBed through pipes
            pipe_t.append("bamToBed " , '--')
            alignment_file_handle = pipe_t.open(None, 'r')
            
    else:
        # There is an actual alignment file
        # that needs to be opened properly
        if alignment_format == 'bam':
            pipe_t.prepend("bamToBed -i {}".format(alignment_file), ".-")
            alignment_file_handle = pipe_t.open(None, 'r')
        else:
            alignment_file_handle = flex_open(alignment_file, 'r')


    return alignment_file_handle


#####################################################################

def flex_open(file, mode = "rt"):
    file_type = type(file)

    if issubclass(file_type, IOBase):
        return file
    elif file_type == str:
        return open_by_extension(file, mode)

#####################################################################

def flex_out(file, mode = "wt"):
    if not file:
        return sys.stdout
    else:
        return flex_open(file , mode)
