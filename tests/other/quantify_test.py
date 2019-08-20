# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import numpy as np
import h5py


from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *
from ribopy.core.quantify import quantify_experiment


import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)
from test_data import *
###########################################

TRANSCRIPT_LENGTHS=\
"""GAPDH    10
VEGFA   5
FLT 22
MYC 1"""


TRANSCRIPT_ANNOTATION=\
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
###########################################

TRANSCRIPT_LENGTHS_2=\
"""GAPDH    90
VEGFA   12
FLT 40
MYC 1
P53 85"""

READ_SET_2=\
"""GAPDH 2  4  read_UTR5_1 0   +
GAPDH 3  5  read_UTR5_2 0   +
GAPDH 10  14  read_UTR5_3 0   +
GAPDH 10  14  read_UTR5_4 0   +
GAPDH 11  14  read_UTR5_5 0   +
GAPDH 12  14  read_UTR5_6 0   +
GAPDH 14  25  read_UTR5_7 0   +
GAPDH 15  20  read_UTR5_junc_1 0   +
GAPDH 16  20  read_UTR5_junc_2 0   +
GAPDH 17  20  read_UTR5_junc_3 0   +
GAPDH 20  25  read_UTR5_junc_4 0   +
GAPDH 23  25  read_UTR5_junc_5 0   +
GAPDH 24  30  read_CDS_1 0   +
GAPDH 30  35  read_CDS_2 0   +
GAPDH 44  46  read_CDS_3 0   +
GAPDH 45  55  read_UTR3_junc_1 0   +
GAPDH 53  57  read_UTR3_junc_2 0   +
GAPDH 54  57  read_UTR3_1 0   +
VEGFA 0  1  read_UTR5_junc_1 0  +
VEGFA 3  5  read_UTR5_junc_2 0  +
VEGFA 4  11  read_CDS_1 0  +
VEGFA 4  5  read_CDS_2 0  +
VEGFA 11  12  read_UTR3_junc_1 0  +
FLT 36  40  read_UTR3_junc_1 0 +
MYC 0  1  read_UTR5_junc_1 0  +
P53 5  15  read_UTR5_1 0  +
P53 25  35  read_CDS_1 0  +"""

ANNOTATION_2 = [\
               ([0, 20], [20, 50], [50, 90]),\
               ([0, 0], [0, 10], [10, 12]),\
               ([0, 15], [15, 40], [40, 40]),\
               ([0, 0], [0, 1], [1, 1]),\
               ([0, 20], [20, 55], [55, 85])
              ]

GAPDH_counts = ( 7, 5, 3, 2, 1 )
VEGFA_counts = ( 0, 2, 2, 1, 0 )
FLT_counts   = ( 0, 0, 0, 1, 0 )
MYC_counts   = ( 0, 1, 0, 0, 0 )
P53_counts   = ( 1, 0, 1, 0, 0 )

expected_counts = ( GAPDH_counts, VEGFA_counts, FLT_counts,
                    MYC_counts, P53_counts )

REF_LEN_FILE_2    = "ref_len_2.txt"
ANNOTATION_FILE_2 = "annotation_2.bed"
ALIGNMENT_FILE_2  = "alignment_2.bed"

###########################################
###########################################


LEFT_SPAN  = 5
RIGHT_SPAN = 3

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
    ref_lengths = tuple( map( int, ref_lengths ) )
    return (ref_names, ref_lengths)


class TestCreate(unittest.TestCase):

    def setUp(self):
        self.tmp_files         = list()
        self.ref_len_file      = StringIO(TRANSCRIPT_LENGTHS)
        self.alignment_file_1  = StringIO(READ_SET_1)
        self.ref_len_file_2    = StringIO(TRANSCRIPT_LENGTHS_2)
        self.alignment_file_2  = StringIO(READ_SET_2)


    def tearDown(self):
        pass


    def test_quantify_metagene(self):
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS)

        ribo_profile = quantify_experiment( 
                          input_reads_file   = self.alignment_file_1, 
                          ref_names          = ref_names, 
                          ref_lengths        = ref_lengths, 
                          region_coordinates = ANNOTATION, 
                          metagene_radius    = 2, 
                          left_span          = LEFT_SPAN, 
                          right_span         = RIGHT_SPAN )

        start_site_coverage = ribo_profile[REF_DG_START_SITE_COV]
        comparison          = (ACTUAL_START_SITE_COVERAGE == \
                               start_site_coverage)
        self.assertTrue( np.all( comparison ) )

        stop_site_coverage = ribo_profile[REF_DG_STOP_SITE_COV]
        comparison         = (ACTUAL_STOP_SITE_COVERAGE ==\
                              stop_site_coverage)
        self.assertTrue( np.all( comparison ) )

    def test_quantify_regions(self):
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS_2)

        ribo_profile = quantify_experiment( 
                          input_reads_file   = self.alignment_file_2, 
                          ref_names          = ref_names, 
                          ref_lengths        = ref_lengths, 
                          region_coordinates = ANNOTATION_2, 
                          metagene_radius    = 2, 
                          left_span          = LEFT_SPAN, 
                          right_span         = RIGHT_SPAN )

        rg = ribo_profile[REF_DG_REGION_COUNTS]

        for i in range(5):
            comparison = (rg[i,: ] == expected_counts[i] )
            self.assertTrue(np.all(comparison) ) 
  

if __name__ == '__main__':
        
    unittest.main()
