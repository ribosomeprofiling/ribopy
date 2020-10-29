# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO, BytesIO

import pandas as pd
import h5py


from ribopy import create
from ribopy.core.get_gadgets import *
from ribopy.settings import *

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)
from test_data import *

###########################################


class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_len_file    = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file  = StringIO(READ_SET_1)
        self.handle          = h5py.File(BytesIO(), "w")

        create.create_ribo(
            self.handle,
            experiment_name    = "merzifon",
            alignment_file  = self.alignment_file,
            metagene_radius = 3,
            length_min      = 2,
            length_max      = 5,
            left_span       = 3,
            right_span      = 2,
            reference_name  = "appris_human_v2",
            lengths_file    = self.ref_len_file,
            annotation_file = self.annotation_file)

    def tearDown(self):
        self.handle.close()


    def test_get_reference_names(self):
        ref_names = get_reference_names(self.handle)
        self.assertListEqual( list(ref_names), ["GAPDH","VEGFA","MYC","BRCA"])

    def test_get_reference_lengths(self):
        ref_lengths = get_reference_lengths(self.handle)
        self.assertListEqual( list(ref_lengths), [1290, 875, 1462, 565])

    def test_get_number_of_references(self):
        number_of_refs = get_number_of_references(self.handle)
        self.assertEqual(number_of_refs, 4)

    def test_get_region_boundaries(self):
        all_boundaries = get_region_boundaries(self.handle)
        self.assertEqual( all_boundaries[1][0][1], 0 )
        self.assertEqual( all_boundaries[3][2][0], 561 )

    def test_get_experiment_names(self):
        self.handle[EXPERIMENTS_name].create_group("Antalya")
        libnames = get_experiment_names(self.handle)
        self.assertTrue( "Antalya" in libnames)
        self.assertTrue( "merzifon" in libnames)

    def test_get_read_length_range(self):
        len_lower , len_upper = get_read_length_range( self.handle )
        self.assertEqual(len_lower, 2)
        self.assertEqual(len_upper, 5)



if __name__ == '__main__':

    unittest.main()
