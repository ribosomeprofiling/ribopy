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

class TestCreate(ApiTestBase):
        
    def test_metagene_empty_exp_list(self):
        
        metagene_start = self.sample_ribo.get_metagene(
                                               site_type      = "start",
                                               sum_lengths    = True,
                                               sum_references = True )
                                               
        self.assertTrue( np.allclose(metagene_start.loc["merzifon"], 
                         [7, 17, 12, 5, 2] )  )
                         
        self.assertTrue( np.allclose(metagene_start.loc["ankara"], 
                         [0, 0, 0, 0, 0] )  )
                         
        metagene_stop = self.sample_ribo.get_metagene(
                                               site_type      = "stop",
                                               sum_lengths    = True,
                                               sum_references = True )
                       
        self.assertTrue( np.allclose(metagene_stop.loc["merzifon"], 
                         [10, 6, 23, 2, 0] )  )
                         
        self.assertTrue( np.allclose(metagene_stop.loc["ankara"], 
                         [0, 0, 0, 0, 0] )  )
                         
    def test_metagene_nonexisting_exp(self):
                         
        with self.assertRaises(RiboBaseError) as context:
            metagene_err = self.sample_ribo.get_metagene(
                                   site_type      = "stop",
                                   sum_lengths    = True,
                                   sum_references = True,
                                   experiments    = ["NONEXISTING_EXP"] )
                                   
    def test_metagene_individual_exp(self):
        metagene_start = self.sample_ribo.get_metagene(
                                               site_type      = "start",
                                               sum_lengths    = True,
                                               sum_references = True,
                                               experiments    = "ankara" )
                                               
        self.assertTrue( np.allclose(metagene_start.loc["ankara"], 
                         [0, 0, 0, 0, 0] )  )
                         
        metagene_stop = self.sample_ribo.get_metagene(
                                               site_type      = "stop",
                                               sum_lengths    = True,
                                               sum_references = True,
                                               experiments    = ["ankara"] )
                                               
        self.assertTrue( np.allclose(metagene_stop.loc["ankara"], 
                         [0, 0, 0, 0, 0] )  )
                         
    def test_metagene_exp_list(self):
        metagene_start = self.sample_ribo.get_metagene(
                               site_type      = "start",
                               sum_lengths    = True,
                               sum_references = True,
                               experiments    = ["ankara", "merzifon"] )
                               
        self.assertTrue( np.allclose(metagene_start.loc["merzifon"], 
                         [7, 17, 12, 5, 2] )  )
                         
        self.assertTrue( np.allclose(metagene_start.loc["ankara"], 
                         [0, 0, 0, 0, 0] )  )
                         
    def test_metagene_ind_transcript(self):
        metagene_start = self.sample_ribo.get_metagene(
                               site_type      = "start",
                               sum_lengths    = True,
                               sum_references = False,
                               experiments    = ["merzifon"] )
                               
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "GAPDH"], 
                                    [7, 7, 6, 2, 2]))
                                  
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "VEGFA"], 
                                    [0, 10, 6, 3, 0]))
                                    
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "MYC"], 
                                    [0, 0, 0, 0, 0]))
                                    
    def test_metagene_ind_length(self):
        metagene_start = self.sample_ribo.get_metagene(
                               site_type      = "start",
                               sum_lengths    = False,
                               range_lower    = 2,
                               range_upper    = 2, 
                               sum_references = False,
                               experiments    = ["merzifon"] )
        
                   
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "GAPDH", 2], 
              ACTUAL_START_SITE_COVERAGE_length_2[0]))
              
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "VEGFA", 2], 
              ACTUAL_START_SITE_COVERAGE_length_2[1]))
              
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "MYC", 2], 
              ACTUAL_START_SITE_COVERAGE_length_2[2]))
              
    def test_metagene_ind_length_range(self):
        metagene_start = self.sample_ribo.get_metagene(
                               site_type      = "start",
                               sum_lengths    = False,
                               range_lower    = 2,
                               range_upper    = 4, 
                               sum_references = False,
                               experiments    = ["merzifon"] )
        
                   
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "GAPDH", 3], 
              ACTUAL_START_SITE_COVERAGE_length_3[0]))
              
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "VEGFA", 3], 
              ACTUAL_START_SITE_COVERAGE_length_3[1]))
              
        self.assertTrue(np.allclose(metagene_start.loc["merzifon", "MYC", 3], 
              ACTUAL_START_SITE_COVERAGE_length_3[2]))
        


        
if __name__ == '__main__':
        
    unittest.main()
