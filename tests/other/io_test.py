# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
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

from ribopy.io.separate_by_length import *

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

# Revise the reads !!!!
READ_SET_1=\
"""GAPDH 2  4000  read_UTR5_1 0   +
GAPDH 3  5000  read_UTR5_2 0   +
GAPDH 10  14000  read_UTR5_3 0   +
GAPDH 10  14000  read_UTR5_4 0   +
GAPDH 11  14000  read_UTR5_5 0   +
GAPDH 12  14000  read_UTR5_6 0   +
GAPDH 14  15  read_UTR5_7 0   +
GAPDH 15  20000  read_UTR5_junc_1 0   +
GAPDH 16  18  read_UTR5_junc_2 0   +
GAPDH 17  20000  read_UTR5_junc_3 0   +
GAPDH 20  22  read_UTR5_junc_4 0   +
GAPDH 23  24  read_UTR5_junc_5 0   +
GAPDH 24  30000  read_CDS_1 0   +
GAPDH 30  32  read_CDS_2 0   +
GAPDH 44  46  read_CDS_3 0   +
GAPDH 45  55000  read_UTR3_junc_1 0   +
GAPDH 53  57000  read_UTR3_junc_2 0   +
GAPDH 54  57000  read_UTR3_1 0   +
VEGFA 0  1000  read_UTR5_junc_1 0  +
VEGFA 3  4  read_UTR5_junc_2 0  +
VEGFA 4  7  read_CDS_1 0  +
VEGFA 4  7  read_CDS_2 0  +
VEGFA 11  12000  read_UTR3_junc_1 0  +
FLT 36  39  read_UTR3_junc_1 0 +
MYC 0  1000  read_UTR5_junc_1 0  +
P53 5  8  read_UTR5_1 0  +
P53 25  27  read_CDS_1 0  +"""

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

### Expected Region Counts (for all read lengths)
# Order is UTR5, UTR5_junc, CDS, UTR3_junc, UTR3

GAPDH_counts_all = ( 7, 5, 3, 2, 1 )
VEGFA_counts_all = ( 0, 2, 2, 1, 0 )
FLT_counts_all   = ( 0, 0, 0, 1, 0 )
MYC_counts_all   = ( 0, 1, 0, 0, 0 )
P53_counts_all   = ( 1, 0, 1, 0, 0 )

expected_counts_all_lengths = ( GAPDH_counts_all, 
                                VEGFA_counts_all, FLT_counts_all,
                                MYC_counts_all, P53_counts_all )

GAPDH_counts_1_to_3 = ( 1, 3, 2, 0, 0 )
VEGFA_counts_1_to_3 = ( 0, 1, 2, 0, 0 )
FLT_counts_1_to_3   = ( 0, 0, 0, 1, 0 )
MYC_counts_1_to_3   = ( 0, 0, 0, 0, 0 )
P53_counts_1_to_3   = ( 1, 0, 1, 0, 0 )

expected_counts_1_to_3 = ( GAPDH_counts_1_to_3, 
                           VEGFA_counts_1_to_3, FLT_counts_1_to_3,
                           MYC_counts_1_to_3, P53_counts_1_to_3 )

GAPDH_counts_2_to_3 = ( 0, 2, 2, 0, 0 )
VEGFA_counts_2_to_3 = ( 0, 0, 2, 0, 0 )
FLT_counts_2_to_3   = ( 0, 0, 0, 1, 0 )
MYC_counts_2_to_3   = ( 0, 0, 0, 0, 0 )
P53_counts_2_to_3   = ( 1, 0, 1, 0, 0 )

expected_counts_2_to_3 = ( GAPDH_counts_2_to_3, 
                           VEGFA_counts_2_to_3, FLT_counts_2_to_3,
                           MYC_counts_2_to_3, P53_counts_2_to_3 )

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


class TestSeparate(unittest.TestCase):

    def setUp(self):
        self.handle = h5py.File(BytesIO())
        self.files_to_be_removed = list()

    def tearDown(self):
        self.handle.close()
        [os.remove(f) for f in self.files_to_be_removed]

    def test_separate_by_length_1_3(self):
        input_stream = StringIO(READ_SET_1)

        separated_reads = separate_by_length(input_stream, length_min = 1, length_max = 3)
        self.files_to_be_removed += separated_reads

        self.assertEqual( len(separated_reads), 3 )

        with open(separated_reads[1], "r") as sep_stream:
            length_2_lines = sep_stream.readlines()
            self.assertEqual(len(length_2_lines), 5)
            parsed_lines = [ k.strip().split() for k in length_2_lines ]
            self.assertListEqual( parsed_lines[1] , 
                     ['GAPDH', '20', '22', 'read_UTR5_junc_4', '0', '+'] )

        with open(separated_reads[0], "r") as sep_stream:
            length_1_lines = sep_stream.readlines()
            self.assertEqual(len(length_1_lines), 3)
            parsed_lines = [ k.strip().split() for k in length_1_lines ]
            self.assertListEqual( parsed_lines[0] ,
                     ['GAPDH', '14', '15', 'read_UTR5_7', '0', '+'])    


if __name__ == '__main__':
        
    unittest.main()
