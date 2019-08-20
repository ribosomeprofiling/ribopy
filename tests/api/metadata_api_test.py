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

class TestHasMetadata_1(ApiTestBase):

    def test_has_metadata(self):
        self.assertTrue(not self.sample_ribo.has_metadata("ankara") )
        self.assertTrue(self.sample_ribo.has_metadata("merzifon") )
        self.assertTrue(self.sample_ribo.has_metadata() )
        
    def test_has_metadata_nonexisting_exp(self):
        with self.assertRaises(ExperimentDoesntExist) as exception_context:
            self.assertTrue(not self.sample_ribo.has_metadata("nonexisting") )
        
class TestGetMetadata_1(ApiTestBase):
    
    def test_get_ribo_metadata(self):
        ribo_metadata = self.sample_ribo.get_metadata()
        self.assertEqual( ribo_metadata.get("aligner"), "bowtie2" )
        self.assertEqual( ribo_metadata.get("deduplicated"), True )
        self.assertEqual( ribo_metadata.get("mapq_threshold"), 3 )
        
    def test_get_exp_metadata(self):
        exp_metadata = self.sample_ribo.get_metadata("merzifon")
        self.assertEqual( exp_metadata.get("cell_line") , "HeLa" )
        self.assertEqual( exp_metadata.get("digestion_enzyme") , "HindIII" )
        self.assertEqual( exp_metadata.get("digestion_duration") , "5 min" )
        self.assertEqual( exp_metadata.get("link") , 
                         "https://www.encodeproject.org/" )
                         
    def test_get_nonexistent_exp_metadata(self):
        with self.assertRaises(RiboBaseError) as exception_context:
            exp_metadata = self.sample_ribo.get_metadata("nonexistent-exp")
        
if __name__ == '__main__':
        
    unittest.main()
