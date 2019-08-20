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
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *
from ribopy.core.aggregate_by_length import *

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)

from multilength_test_data import *
from ribopy.core import create_experiment

###########################################

NPROCESS = 4

class TestCreate(unittest.TestCase):
    def setUp(self):

        self.len_file        = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file  = StringIO(READ_SET_1)

        self.handle = h5py.File(BytesIO())

        create.create_ribo(self.handle, 
                experiment_name        = "merzifon", 
                alignment_file      = self.alignment_file,
                reference_name      = "hg38",
                lengths_file        = self.len_file, 
                annotation_file     = self.annotation_file,
                metagene_radius     = METAGENE_RADIUS, 
                left_span = LEFT_SPAN, right_span = RIGHT_SPAN,
                length_min = 2, length_max = 5, 
                tmp_file_prefix     = "")
        
    def tearDown(self):
        self.handle.close()


    ### T E S T    S U M S  ################################    

    def test_aggregate_region_sum(self):
        aggregate_2_3_region = aggregate_region_counts(
                                 self.handle, range_lower = 2, range_upper =3,
                                 sum_values = True )
        expected_reg_count_len_2_3 = expected_counts_length_2 + \
                                       expected_counts_length_3

        # print(aggregate_2_3_region, "\n-----\n" , expected_reg_count_len_2_3 )
        self.assertTrue( np.all(aggregate_2_3_region ==\
                                expected_reg_count_len_2_3) )

        expected_reg_count_len_2_3_4_5 = expected_counts_length_2 + \
                                       expected_counts_length_3 + \
                                       expected_counts_length_4 + \
                                       expected_counts_length_5

        aggregate_2_5_region = aggregate_region_counts(
                                 self.handle, range_lower = 2, range_upper =5,
                                 sum_values = True )

        self.assertTrue( np.all(expected_reg_count_len_2_3_4_5 ==\
                                aggregate_2_5_region) )

        aggregate_2_region = aggregate_region_counts(
                                 self.handle, range_lower = 2, range_upper =2,
                                 sum_values      = True,
                                 experiment_list = ["merzifon", "merzifon"] )
        self.assertTrue( np.all(aggregate_2_region ==\
                                expected_counts_length_2) )


    def test_aggregate_start_site_sum(self):
        aggregate_2_3_start = aggregate_start_site(
                                 self.handle, range_lower = 2, range_upper =3,
                                 sum_values = True )
        expected_2_3_start = ACTUAL_START_SITE_COVERAGE_length_2 + \
                             ACTUAL_START_SITE_COVERAGE_length_3

        self.assertTrue( np.all(aggregate_2_3_start == expected_2_3_start)  )

        aggregate_2_start = aggregate_start_site(
                                 self.handle, range_lower = 2, range_upper =2,
                                 sum_values = True )

        self.assertTrue( np.all(aggregate_2_start ==\
                                ACTUAL_START_SITE_COVERAGE_length_2)  )

        aggregate_2_to_5_start = aggregate_start_site(
                                 self.handle, range_lower = 2, range_upper =5,
                                 sum_values = True )

        expected_2_to_5_start = ACTUAL_START_SITE_COVERAGE_length_2 + \
                             ACTUAL_START_SITE_COVERAGE_length_3 + \
                             ACTUAL_START_SITE_COVERAGE_length_4 + \
                             ACTUAL_START_SITE_COVERAGE_length_5

        self.assertTrue( np.all(aggregate_2_to_5_start == \
                                expected_2_to_5_start)  )


    def test_aggregate_stop_site_sum(self):
        aggregate_2_3_stop = aggregate_stop_site(
                                 self.handle, range_lower = 2, range_upper =3,
                                 sum_values = True )

        expected_2_3_stop = ACTUAL_STOP_SITE_COVERAGE_length_2 + \
                             ACTUAL_STOP_SITE_COVERAGE_length_3

        self.assertTrue( np.all(aggregate_2_3_stop == expected_2_3_stop)  )


    ### Test GROUPINGS #######################
    def test_aggregate_region_group(self):
        aggregate_2_3_region = aggregate_region_counts(
                                 self.handle, range_lower = 2, range_upper =3,
                                 sum_values = False )

        length_col  = np.repeat( range(2,3+1), len( expected_counts_length_2 ) )
        length_col =  np.array( length_col, ndmin=2 )

        combined_counts = np.concatenate( (expected_counts_length_2,
                                           expected_counts_length_3) )

        expected_array = np.concatenate( (length_col.T, combined_counts), axis = 1 )

        array_comparison = (expected_array == aggregate_2_3_region)
        self.assertTrue( np.all(array_comparison) )

    def test_aggregate_start_site_group(self):
        aggregate_4_5_start = aggregate_start_site(
                                 self.handle, range_lower = 4, range_upper =5,
                                 sum_values = False )

        length_col  = np.repeat( range(4,5+1), 
                                 len( ACTUAL_START_SITE_COVERAGE_length_4 ) )
        length_col =  np.array( length_col, ndmin=2 )

        combined_counts = np.concatenate( (ACTUAL_START_SITE_COVERAGE_length_4,
                                           ACTUAL_START_SITE_COVERAGE_length_5) )

        expected_array = np.concatenate( (length_col.T, combined_counts), axis = 1 )
        array_comparison = (expected_array == aggregate_4_5_start)
        self.assertTrue( np.all(array_comparison) )

    def test_aggregate_stop_site_group(self):
        aggregate_2_5_stop = aggregate_stop_site(
                                 self.handle, range_lower = 2, range_upper =5,
                                 sum_values = False )

        length_col  = np.repeat( range(2,5+1), len( ACTUAL_STOP_SITE_COVERAGE_length_4 ) )
        length_col =  np.array( length_col, ndmin=2 )

        combined_counts = np.concatenate( (ACTUAL_STOP_SITE_COVERAGE_length_2,
                                           ACTUAL_STOP_SITE_COVERAGE_length_3,
                                           ACTUAL_STOP_SITE_COVERAGE_length_4,
                                           ACTUAL_STOP_SITE_COVERAGE_length_5) )

        expected_array = np.concatenate( (length_col.T, combined_counts), axis = 1 )
        array_comparison = (expected_array == aggregate_2_5_stop)
        self.assertTrue( np.all(array_comparison) )



if __name__ == '__main__':
        
    unittest.main()
