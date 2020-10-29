# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO, BytesIO

import numpy as np
import h5py

from test_data import *
from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.metagene import find_site_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)
from test_data import *

###########################################

TRANSCRIPT_LENGTHS=\
"""GAPDH    10
VEGFA   5
FLT 22
MYC 2"""


TRANSCRIPT_ANNOTATION=\
"""GAPDH    0   4  UTR5    0   -
GAPDH   4  7    CDS 0   -
GAPDH   7    10 UTR3   0   -
VEGFA   0   3 CDS 0   +
VEGFA   3 5 UTR3    0   +
FLT   0 6 UTR5    0   +
FLT   6 15 CDS    0   +
FLT   15 22 UTR3    0   +
MYC 0  1    CDS 0   -
MYC 1    2    UTR3    0   -"""

READ_SET_1=\
"""GAPDH 2  4  read_1 0   +
MYC 0  1  read_2 0 +
VEGFA 0  5  read_30 0 +
VEGFA 0  6  read_31 0 +
VEGFA 4  5  read_3 0 +
MYC 0  1  read_4 0 +
MYC 0  1  read_5 0 +
VEGFA 4  5  read_6 0 +
VEGFA 1  2  read_7 0 +
FLT 4  10  read_8 0 +
FLT 4  6  read_9 0 +
FLT 4  8  read_10 0 +
FLT 4  7  read_11 0 +
FLT 4  9  read_12 0 +
FLT 4  5  read_13 0 +
FLT 5  6  read_14 0 +
FLT 5  6  read_15 0 +
FLT 5  7  read_16 0 +
FLT 5  8  read_17 0 +
FLT 6  9  read_18 0 +
FLT 7  8  read_19 0 +
FLT 7  9  read_20 0 +
FLT 7  10  read_21 0 +
FLT 7  11  read_22 0 +
FLT 7  12  read_23 0 +
FLT 8  9  read_24 0 +
FLT 8  11  read_240 0 +
FLT 9  19  read_25 0 +
FLT 10  12  read_26 0 +
FLT 11  15  read_27 0 +
FLT 12  14  read_28 0 +
FLT 13  15  read_29 0 +
FLT 13  16  read_30 0 +
FLT 13  17  read_31 0 +
FLT 13  17  read_32 0 +
FLT 14  17  read_33 0 +
FLT 14  17  read_34 0 +
FLT 16  17  read_35 0 +
FLT 16  17  read_36 0 +
FLT 16  17  read_37 0 +
FLT 17  19  read_38 0 +"""

"""
Expected start site distribution, for radius=2:

GAPDH: [1,0,0,0,0]
VEGFA: [0, 0, 2, 1, 0]
FLT1: [6, 4, 1, 5, 2]
MYC: [0, 0, 3, 0, 0]
"""

"""
Expected stop site distribution, for radius = 2:

GAPDH: [0, 0, 0, 0, 0]
VEGFA: [1, 0, 0, 2, 0]
MYC: [0, 3, 0, 0, 0 ]
FLT1: [4, 2, 0, 3, 1]
"""

ACTUAL_START_SITE_COVERAGE = np.array( [[1, 0, 0, 0, 0],
                                [0, 0, 2, 1, 0],
                                [6, 4, 1, 5, 2],
                                [0, 0, 3, 0, 0]] )

ACTUAL_STOP_SITE_COVERAGE = np.array( [ [0, 0, 0, 0, 0],
                                        [1, 0, 0, 2, 0],
                                        [4, 2, 0, 3, 1],
                                        [0, 3, 0, 0, 0]] )

###########################################
"""GAPDH    0   4  UTR5    0   -
GAPDH   4  7    CDS 0   -
GAPDH   7    10 UTR3   0   -
VEGFA   0   3 CDS 0   +
VEGFA   3 5 UTR3    0   +
FLT   0 6 UTR5    0   +
FLT   6 15 CDS    0   +
FLT   15 20 UTR3    0   +
MYC 0  1    CDS 0   -
MYC 1    1    UTR3    0   -"""


ANNOTATION = [\
               ([0, 4], [4, 7], [7, 10]),\
               ([0, 0], [0, 3], [3, 5]),\
               ([0, 6], [6, 15], [15, 20]),\
               ([0,0], [0, 1], [1,1])
              ]


def _get_transcripts(file_in_string):
    rows        = file_in_string.split("\n")
    pairs       = tuple( map( lambda x: x.split(), rows ) )
    ref_names   = tuple( map( lambda x: x[0], pairs ) )
    ref_lengths = tuple( map( lambda x: x[1], pairs ) )
    ref_lengths = tuple (map( int, ref_lengths ) )
    return (ref_names, ref_lengths)


class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_len_file    = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file  = StringIO(READ_SET_1)
        self.handle          = h5py.File(BytesIO(), "w")

    def tearDown(self):
        self.handle.close()


    def test_start_site_coverage_solo(self):
        input_stream = StringIO(READ_SET_1)
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS)

        coverage = find_coverage(input_stream,
                                 ref_names , ref_lengths )

        site_coverage = find_site_coverage(coverage   = coverage,
                                           radius     = 2,
                                           annotation = ANNOTATION,
                                           site_type  = "start")

        comparison = (ACTUAL_START_SITE_COVERAGE == site_coverage)
        self.assertTrue( np.all( comparison ) )

    def test_stop_site_coverage_solo(self):
        input_stream = StringIO(READ_SET_1)
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS)

        coverage = find_coverage(input_stream,
                                 ref_names , ref_lengths )

        site_coverage = find_site_coverage(coverage   = coverage,
                                           radius     = 2,
                                           annotation = ANNOTATION,
                                           site_type  = "stop")

        comparison = (ACTUAL_STOP_SITE_COVERAGE == site_coverage)
        self.assertTrue( np.all( comparison ) )

    def test_start_site_coverage_from_ribo(self):
        create.initialize(self.handle, "merzifon",
                          reference_name = "appris_human_v2")
        create.set_reference_names_and_lengths(
                                       self.handle,
                                       self.ref_len_file)
        create.set_annotation(
            h5_handle       = self.handle,
            annotation_lines = TRANSCRIPT_ANNOTATION.split("\n") )

        ribo_annotation = get_region_boundaries(self.handle)
        ref_names       = get_reference_names(self.handle)
        ref_lengths     = get_reference_lengths(self.handle)
        input_stream    = StringIO(READ_SET_1)
        coverage        = find_coverage(input_stream,
                                        ref_names , ref_lengths )

        site_coverage = find_site_coverage(coverage   = coverage,
                                           radius     = 2,
                                           annotation = ribo_annotation,
                                           site_type  = "start")

        comparison = (ACTUAL_START_SITE_COVERAGE == site_coverage)
        self.assertTrue( np.all( comparison ) )


if __name__ == '__main__':

    unittest.main()
