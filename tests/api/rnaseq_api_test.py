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
from ribopy.rnaseq import set_rnaseq, get_rnaseq

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from multilength_test_data import *
from api_test_base import ApiTestBase
from rnaseq_data import *

####################################################################

NPROCESS = 1

class TestRNASEQ(ApiTestBase):

    def setUp(self):
        super().setUp()
        # We need to init the ribo object in write mode
        # Because we are going to add tnaseq data to it.
        self.sample_ribo = Ribo( self.merged_io, file_mode = "r+" )
        
        self.rnaseq_file_1 = StringIO( RNASEQ_READS )
        self.rnaseq_file_2 = StringIO( RNASEQ_READS_2 )
        
        set_rnaseq(ribo_handle   = self.sample_ribo._handle, 
                   name          = "merzifon", 
                   rnaseq_reads  = self.rnaseq_file_1, 
                   format        = "bed",
                   rnaseq_counts = None)
        
        # We need to re-initialize the ribo object
        # so that the changes take effect.
        self.sample_ribo._handle.close()
        self.sample_ribo = Ribo( self.merged_io, file_mode = "r+" )
    
    def test_has_rnaseq(self):
        self.assertTrue(self.sample_ribo.has_rnaseq("merzifon")) 
        
        self.assertTrue(not self.sample_ribo.has_rnaseq("ankara"))
        
    def test_rnaseq_info(self):
        self.assertTrue(self.sample_ribo.\
                info["experiments"]["merzifon"]["RNA-Seq"])
        self.assertTrue(not self.sample_ribo.\
               info["experiments"]["ankara"]["RNA-Seq"])  
               
    def test_get_rnaseq_single(self):
        rnaseq_df = self.sample_ribo.get_rnaseq("merzifon")
        
        self.assertTrue(np.allclose( [1, 0, 2, 1, 3], 
                                     rnaseq_df.loc["merzifon", "GAPDH"] ))  
        self.assertTrue(np.allclose( [0, 1, 3, 0, 1], 
                                     rnaseq_df.loc["merzifon", "VEGFA"] )) 
        self.assertTrue(np.allclose( [0, 1, 0, 0, 2], 
                                     rnaseq_df.loc["merzifon", "MYC"] ))
                                     
        self.assertTrue( np.isclose(rnaseq_df.loc["merzifon", "GAPDH"][CDS_name] , 
                                    2 ) )
                                    
    def test_get_rnaseq_from_nornaseq_exp(self):
        with self.assertRaises(NORNASEQ) as exception_context:
            rnaseq_df = self.sample_ribo.get_rnaseq("ankara")
                                    
class TestRNASEQ_2(ApiTestBase):

    def setUp(self):
        super().setUp()
        # We need to init the ribo object in write mode
        # Because we are going to add tnaseq data to it.
        self.sample_ribo = Ribo( self.merged_io, file_mode = "r+" )
        
        self.rnaseq_file_1 = StringIO( RNASEQ_READS )
        self.rnaseq_file_2 = StringIO( RNASEQ_READS_2 )
        
        set_rnaseq(ribo_handle   = self.sample_ribo._handle, 
                   name          = "merzifon", 
                   rnaseq_reads  = self.rnaseq_file_1, 
                   format        = "bed",
                   rnaseq_counts = None)
                   
        set_rnaseq(ribo_handle   = self.sample_ribo._handle, 
                   name          = "ankara", 
                   rnaseq_reads  = self.rnaseq_file_2, 
                   format        = "bed",
                   rnaseq_counts = None)
        
        # We need to re-initialize the ribo object
        # so that the changes take effect.
        self.sample_ribo._handle.close()
        self.sample_ribo = Ribo( self.merged_io, file_mode = "r+" )
    
    def test_has_rnaseq(self):
        self.assertTrue(self.sample_ribo.has_rnaseq("merzifon")) 
        
        self.assertTrue(self.sample_ribo.has_rnaseq("ankara"))
        
               
    def test_get_rnaseq_multiple(self):
        rnaseq_df = self.sample_ribo.get_rnaseq()
        
        self.assertTrue(np.allclose( [1, 0, 2, 1, 3], 
                                     rnaseq_df.loc["merzifon", "GAPDH"] ))  
        self.assertTrue(np.allclose( [0, 1, 3, 0, 1], 
                                     rnaseq_df.loc["merzifon", "VEGFA"] )) 
        self.assertTrue(np.allclose( [0, 1, 0, 0, 2], 
                                     rnaseq_df.loc["merzifon", "MYC"] ))
                                     
        self.assertTrue( np.isclose(rnaseq_df.loc["merzifon", "GAPDH"][CDS_name] , 
                                    2 ) ) 
                                    
        self.assertTrue(np.allclose( [0, 0, 0, 0, 0], 
                                     rnaseq_df.loc["ankara", "GAPDH"] ))  
        self.assertTrue(np.allclose( [0, 0, 0, 2, 0], 
                                     rnaseq_df.loc["ankara", "VEGFA"] )) 
        self.assertTrue(np.allclose( [0, 0, 0, 0, 0], 
                                     rnaseq_df.loc["ankara", "MYC"] ))
                                     
        self.assertTrue( np.isclose(rnaseq_df.loc["ankara", "VEGFA"][UTR3_JUNCTION_name] , 
                                    2 ) )
                                    
                                     

        
if __name__ == '__main__':
        
    unittest.main()
