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
from ribopy.dump import *
from ribopy.core import create_experiment

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)

from multilength_test_data import *


###########################################

NPROCESS = 4

class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_len_file    = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file  = StringIO(READ_SET_1)

        self.handle = h5py.File(BytesIO(), "w")

        create.create_ribo(     self.handle, "merzifon",
                                alignment_file  = self.alignment_file,
                                reference_name  = "appris_human_v2",
                                length_min      = 2,
                                length_max      = 5,
                                metagene_radius = METAGENE_RADIUS,
                                left_span       = LEFT_SPAN,
                                right_span      = RIGHT_SPAN,
                                lengths_file    = self.ref_len_file,
                                store_coverage  = True,
                                annotation_file = self.annotation_file)

    def tearDown(self):
        self.handle.close()

    ### T E S T    S U M S  ################################

    def test_dump_metagene_startsite(self):
        dumped_df = dump_metagene(
                      ribo_handle     = self.handle,
                      site_type       = REF_DG_START_SITE_COV,
                      sum_lengths     = True,
                      sum_references  = False,
                      range_lower     = 2,
                      range_upper     = 4,
                      experiment_list = ["merzifon"] )["merzifon"]

        self.assertTrue( list(dumped_df.iloc[1]) == [0, 9, 6, 3, 0])

    def test_dump_metagene_stopsite(self):
        dumped_df = dump_metagene(
                      ribo_handle     = self.handle,
                      site_type       = REF_DG_STOP_SITE_COV,
                      sum_lengths     = True,
                      sum_references  = False,
                      range_lower     = 2, range_upper = 4,
                      experiment_list = ["merzifon"] )["merzifon"]

        # Start site is different from stop site
        self.assertTrue( not list(dumped_df.iloc[1]) == [0, 9, 6, 3, 0])

        self.assertTrue( list(dumped_df.iloc[2]) == [0, 0, 6, 0, 0])


    def test_dump_metagene_start_all_lengths(self):
        dumped_df = dump_metagene(
                      ribo_handle     = self.handle,
                      site_type       = REF_DG_START_SITE_COV,
                      sum_lengths     = False,
                      sum_references  = False,
                      range_lower     = 2,
                      range_upper     = 3,
                      experiment_list = ["merzifon"] )["merzifon"]

        # Note that the first 3 is for the read length
        self.assertTrue( list(dumped_df.iloc[4] ) == [3,0,3,2,1,0] )


    def test_dump_metagene_start_all_lengths(self):
        dumped_df = dump_region(ribo_handle     = self.handle,
                                region_name     = CDS_name,
                                sum_lengths     = True,
                                sum_references  = False,
                                range_lower     = 2 ,
                                range_upper     = 4,
                                experiment_list = [])["merzifon"]

        self.assertTrue( list(dumped_df[...]) ==  [11, 6 ,3])

    def test_dump_coverage(self):
        len_4_cov = dump_coverage( ribo_handle     = self.handle,
                                   experiment_name = "merzifon",
                                   range_lower     = 4,
                                   range_upper     = 4)


        self.assertEqual(len_4_cov["GAPDH"][7], 1)


        len_3_and_4_cov = dump_coverage(
                                   ribo_handle     = self.handle,
                                   experiment_name = "merzifon",
                                   range_lower     = 3,
                                   range_upper     = 4)
        self.assertEqual(len_3_and_4_cov["VEGFA"][4], 5)

if __name__ == '__main__':

    unittest.main()
