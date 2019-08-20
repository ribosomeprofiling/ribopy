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
        cds_counts_sum_all = \
           self.sample_ribo.get_region_counts( region_name = "CDS")
        
        self.assertEqual( cds_counts_sum_all["merzifon"][0], 26  )
        self.assertEqual( cds_counts_sum_all["ankara"][0], 8  )
        
        cds_counts_sum_ref = \
           self.sample_ribo.get_region_counts( region_name    = "CDS",
                                               sum_lengths    = False)
                                               
        self.assertTrue( np.allclose( cds_counts_sum_ref["merzifon"][...] , 
                                     [3, 6, 11, 6] ) )
                                     
        self.assertTrue( np.allclose( cds_counts_sum_ref["ankara"][...] , 
                                     [0, 8, 0, 0] ) )
        
        cds_counts_sum_len = \
           self.sample_ribo.get_region_counts( region_name    = "CDS",
                                               sum_references = False)
        self.assertTrue( np.isclose( cds_counts_sum_len.loc["MYC"]["merzifon"],
                                   3))
        self.assertTrue( np.isclose( cds_counts_sum_len.loc["GAPDH"]["merzifon"],
                                   13))
        self.assertTrue( np.isclose( cds_counts_sum_len.loc["VEGFA"]["merzifon"],
                                   10))                                          
        
        cds_counts_no_sum = \
           self.sample_ribo.get_region_counts( region_name    = "CDS",
                                               sum_references = False,
                                               sum_lengths    = False)
                                               
        expected_merzifon_no_sum = [2, 1, 0, 3, 2, 1, 6, 3, 2, 2, 4, 0]
        
        self.assertTrue( np.allclose( cds_counts_no_sum["merzifon"][...],
                                      expected_merzifon_no_sum ))
                                      
    
    def test_metagene_specified_exp_list(self):
        utr5_counts_sum_all = \
           self.sample_ribo.get_region_counts( region_name = "UTR5",
                                               experiments = "merzifon")
                                               
        self.assertTrue( "merzifon" in utr5_counts_sum_all.columns)
        self.assertTrue( "ankara" not in utr5_counts_sum_all.columns)
        
        self.assertTrue( np.isclose(utr5_counts_sum_all["merzifon"][0], 4 ) )
        
        utr3_counts_sum_all = \
           self.sample_ribo.get_region_counts( region_name = "UTR3",
                                               experiments = ["merzifon"])
                                               
        self.assertTrue( np.isclose(utr3_counts_sum_all["merzifon"][0], 2 ) )
        
    ################################################################    
    def test_junction_regions(self):
        utr5_junction_sum_all = \
            self.sample_ribo.get_region_counts( region_name = "UTR5_junction",
                                                experiments = "merzifon")
        
        self.assertTrue( np.isclose(utr5_junction_sum_all["merzifon"][0], 45 ) )
        
        utr3_junction_sum_all = \
            self.sample_ribo.get_region_counts( region_name = "UTR3_junction",
                                                experiments = "merzifon")
        
        self.assertTrue( np.isclose(utr3_junction_sum_all["merzifon"][0], 41 ) )
                                                

        
if __name__ == '__main__':
        
    unittest.main()
