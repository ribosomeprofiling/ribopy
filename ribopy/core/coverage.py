# -*- coding: utf-8 -*-
from collections import OrderedDict

import h5py
import numpy as np

from ..settings import *

###############################################################

REF_INDEX        = 0
FIVE_P_END_INDEX = 1


def find_coverage( input_bed_stream,  
                   ref_names , ref_lengths):
    """
    Computes the coverage per transcript.

    Note that the order in the ref_names and the ref_lengths must be compatible.

    This order is also kept in the order of the keys of the
    returned value ( coverage_dict )

    Parameters
    ----------
    input_bed_stream: IO_stream
       Alignemt reads are in bed format.
       This is a handle to the bed file

    ref_names: list(str)
       List of transcript names

    ref_lengths: list(int)
       List of transcript lengths.
       They are in the same order with the ref names.

    Returns
    -------
    coverage_dict: OrderedDict(numpy.array)
       The keys are transcript names in the same
       order as ref_names
       The values are numpy arrays corresponding
       to coverage at the given index position.

    """
    assert len(ref_names) == len(ref_lengths)

    coverage_dict = OrderedDict()

    for ref_name, ref_length in zip(ref_names, ref_lengths):
        coverage_dict[ref_name] = np.zeros( ref_length, 
                                            dtype = TRANSCRIPT_COVERAGE_DT) 

    for entry in input_bed_stream:
        contents = entry.strip().split()
        if len(contents) < 6:
            continue
        this_ref    = contents[REF_INDEX]
        this_five_p = int(contents[FIVE_P_END_INDEX])
        coverage_dict[this_ref][ this_five_p ] += 1

    return coverage_dict