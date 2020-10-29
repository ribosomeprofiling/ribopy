# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import h5py

from ribopy import Ribo
from ribopy import create
from ribopy.merge import merge_ribos
from ribopy.settings import *
from ribopy.core.exceptions import *

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from multilength_test_data import *
from api_test_base import ApiTestBase

####################################################################

class TestHasCoverage(ApiTestBase):

    def test_has_coverage(self):
        self.assertTrue(not self.sample_ribo.has_coverage("merzifon") )
        self.assertTrue( self.sample_ribo.has_coverage("ankara") )

    def test_has_coverage_nonexising_exp(self):
        with self.assertRaises(ExperimentDoesntExist) as exception_context:
            has_coverage = self.sample_ribo.has_coverage("nonexisting")

class TestGetCoverage(ApiTestBase):
    def setUp(self):
        super().setUp()
        self.file_handle_3 =  BytesIO()
        self.ribo_handle_3 = h5py.File( self.file_handle_3, "w" )

        create.create_ribo(
                ribo            = self.ribo_handle_3 ,
                experiment_name = "adana",
                alignment_file  = self.alignment_file_1,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file,
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN,
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                store_coverage  = True,
                nprocess        = 1,
                tmp_file_prefix = "")


        self.ribo_handle_3.close()
        self.adana_ribo = Ribo(self.file_handle_3)

    def test_get_coverage(self):
        coverage = self.adana_ribo.get_coverage("adana")
        self.assertEqual( coverage["GAPDH"][4], 7 )
        self.assertEqual( coverage["VEGFA"][15], 0 )
        self.assertEqual( coverage["MYC"][10], 13 )
        self.assertEqual( coverage["MYC"][0], 2 )

    def test_get_coverage_nonexising_exp(self):
        with self.assertRaises(ExperimentDoesntExist) as exception_context:
            has_coverage = self.sample_ribo.get_coverage("nonexisting")

class TestGetTranscriptCoverage(TestGetCoverage):

    def test_get_transcript_coverage(self):
        """
        all_df = coverage = self.sample_ribo.get_coverage("ankara")
        print(all_df["GAPDH"])

        print("------------------")
        coverage = self.sample_ribo.get_transcript_coverage(
                          experiment = "ankara",
                          transcript = "GAPDH",
                          range_upper = 5,
                          range_lower = 2,
                          sum_lengths = False)
        print(coverage)
        """
        coverage = self.sample_ribo.get_transcript_coverage(
                          experiment = "ankara",
                          transcript = "GAPDH",
                          range_upper = 3,
                          range_lower = 3,
                          sum_lengths = False)


        self.assertEqual(coverage.loc[3][8], 5)

        coverage = self.adana_ribo.get_transcript_coverage(
                          experiment = "adana",
                          transcript = "GAPDH",
                          range_upper = 5,
                          range_lower = 2,
                          sum_lengths = False)

        self.assertTrue( np.all(np.isclose(coverage[4], [2,2,3, 0]) ) )

if __name__ == '__main__':

    unittest.main()
