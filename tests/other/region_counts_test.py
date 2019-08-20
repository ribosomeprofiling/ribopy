# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import numpy as np
import h5py


from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.region_counts import get_extended_boundaries,\
                                           find_region_counts

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
"""GAPDH    90
VEGFA   12
FLT 40
MYC 1
P53 85"""


TRANSCRIPT_ANNOTATION=\
"""GAPDH    0   20  UTR5    0   -
GAPDH   20  50    CDS 0   -
GAPDH   50    90 UTR3   0   -
VEGFA   0   10 CDS 0   +
VEGFA   10 12 UTR3    0   +
FLT   0 15 UTR5    0   +
FLT   15 40 CDS    0   +
MYC 0  1    CDS 0   -
MYC 1    1    UTR3    0   -
P53 0   20  UTR5    0   +
P53 20  55  CDS 0   +
P53 55  85  UTR3    0   +"""


READ_SET_1=\
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

#########################################
#####  Expected EXTENDED ANNOTATION 

GAPDH_regions = ( (0, 15), (15, 24), (24, 45), (45, 54), (54, 90) )
VEGFA_regions = ( (0, 0), (0, 4), (4, 5), (5, 12), (12, 12) )
FLT_regions   = ( (0, 10), (10, 19), (19, 35), (35, 40), (40, 40) )
MYC_regions   = ( (0, 0), (0, 1), (1, 1), (1, 1), (1, 1) )
P53_regions   = ( (0, 15), (15, 24), (24, 50), (50, 59), (59, 85) ) 
extended_boundary_regions = ( GAPDH_regions, VEGFA_regions,
                              FLT_regions, MYC_regions, P53_regions )
#########################################

#########################################

### Expected Region Counts
# Order is UTR5, UTR5_junc, CDS, UTR3_junc, UTR3

GAPDH_counts = ( 7, 5, 3, 2, 1 )
VEGFA_counts = ( 0, 2, 2, 1, 0 )
FLT_counts = ( 0, 0, 0, 1, 0 )
MYC_counts = ( 0, 1, 0, 0, 0 )
P53_counts = ( 1, 0, 1, 0, 0 )

expected_counts = ( GAPDH_counts, VEGFA_counts, FLT_counts,
                    MYC_counts, P53_counts )

#########################################


LEFT_SPAN  = 5
RIGHT_SPAN = 3

ANNOTATION = [\
               ([0, 20], [20, 50], [50, 90]),\
               ([0, 0], [0, 10], [10, 12]),\
               ([0, 15], [15, 40], [40, 40]),\
               ([0, 0], [0, 1], [1, 1]),\
               ([0, 20], [20, 55], [55, 85])
              ]


def _get_transcripts( file_in_string ):
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
        self.handle          = h5py.File(BytesIO())
        

    def tearDown(self):
        self.handle.close()


    def test_get_extended_boundaries(self):
        input_stream = StringIO( READ_SET_1 )
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS)

        coverage = find_coverage(input_stream,  ref_names , ref_lengths )
        extended_boundaries = \
              get_extended_boundaries( annotation = ANNOTATION, 
                                       left_span  = LEFT_SPAN, 
                                       right_span = RIGHT_SPAN )

        for computed, expected in \
              zip(extended_boundaries , extended_boundary_regions):
            self.assertTrue( computed == expected ) 

  
    def test_find_region_counts(self):
        input_stream = StringIO(READ_SET_1)
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS)

        coverage = find_coverage(input_stream,  ref_names , ref_lengths )
        rg = find_region_counts( coverage   = coverage, 
                                 annotation = ANNOTATION, 
                                 left_span  = LEFT_SPAN, 
                                 right_span = RIGHT_SPAN)

        for i in range(5):
            comparison = (rg[i,: ] == expected_counts[i] )
            self.assertTrue(np.all(comparison) ) 


if __name__ == '__main__':
        
    unittest.main()
