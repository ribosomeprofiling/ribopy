# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import pandas as pd
import h5py

from test_data import *
from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths
from ribopy.settings import *


###########################################

TRANSCRIPT_LENGTHS=\
"""GAPDH    10
VEGFA   5
MYC 2"""


TRANSCRIPT_ANNOTATION=\
"""GAPDH    0   4  UTR5    0   -
GAPDH   4  7    CDS 0   -
GAPDH   7    10 UTR3   0   -
VEGFA   0   3 CDS 0   +
VEGFA   3 5 UTR3    0   +
MYC 0  1    CDS 0   -
MYC 1    2    UTR3    0   -"""

READ_SET_1=\
"""GAPDH 2  4  read_1 0   +
MYC 0  1  read_2 0 +
VEGFA 4  5  read_3 0 +
MYC 0  1  read_4 0 +
MYC 0  1  read_5 0 +
VEGFA 4  5  read_6 0 +
VEGFA 1  2  read_7 0 +"""

###########################################


def check_READ_SET_1(self, result):
    self.assertEqual( result["MYC"][0], 3 )
    self.assertEqual( result["VEGFA"][4], 2 )
    self.assertEqual( result["VEGFA"][1], 1 )
    self.assertNotEqual( result["VEGFA"][0], 2 )
    self.assertEqual( result["GAPDH"][2], 1 )    

def _get_transcripts(file_in_string):
    rows = file_in_string.split("\n")
    pairs = tuple( map( lambda x: x.split(), rows ) )
    ref_names = tuple( map( lambda x: x[0], pairs ) )
    ref_lengths = tuple( map( lambda x: x[1], pairs ) )
    ref_lengths = tuple (map( int, ref_lengths ) )
    return (ref_names, ref_lengths)


class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_len_file     = StringIO( TRANSCRIPT_LENGTHS )
        self.annotation_file  = StringIO( TRANSCRIPT_ANNOTATION )
        self.alignment_file   = StringIO(READ_SET_1)

        self.handle = h5py.File(BytesIO())

    def tearDown(self):
        self.handle.close()

    def test_coverage_solo(self):
        input_stream = StringIO(READ_SET_1)
        ref_names , ref_lengths = _get_transcripts(TRANSCRIPT_LENGTHS)
        result = find_coverage( input_stream, ref_names , ref_lengths )
        check_READ_SET_1(self, result)

    def test_coverage_from_ribo(self):
        create.initialize(self.handle, "merzifon", 
                          reference_name = "appris_human_v2")
        create.set_reference_names_and_lengths(self.handle, self.ref_len_file)
        create.set_annotation( h5_handle        = self.handle, 
                               annotation_lines = TRANSCRIPT_ANNOTATION.split("\n") )

        input_stream = StringIO(READ_SET_1)
        ref_names    = get_reference_names(self.handle)
        ref_lengths  = get_reference_lengths(self.handle)
        result       = find_coverage( input_stream, ref_names , ref_lengths )
        check_READ_SET_1(self, result)
        

if __name__ == '__main__':
    unittest.main()
