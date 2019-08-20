# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import h5py

from ribopy import Ribo
from ribopy import create
from ribopy.merge import merge_ribos
from ribopy.settings import *

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from multilength_test_data import *
from api_test_base import ApiTestBase

####################################################################

class TestCreate(ApiTestBase):

    def test_ribo_info(self):
        
        ribo_info   = self.sample_ribo.info
        
        experiments =  tuple( ribo_info[EXPERIMENTS_name].keys() )
        self.assertTrue("ankara" in experiments)
        self.assertTrue("merzifon" in experiments)
        self.assertEqual( len(experiments), 2 )
        
        self.assertEqual( ribo_info[EXPERIMENTS_name]["ankara"]["Coverage"], True )
        self.assertEqual( ribo_info[EXPERIMENTS_name]["merzifon"]["Coverage"], False )
        
        self.assertEqual( ribo_info[EXPERIMENTS_name]["ankara"]["RNA-Seq"], False )
        self.assertEqual( ribo_info[EXPERIMENTS_name]["merzifon"]["RNA-Seq"], False )
        
        self.assertEqual( ribo_info[EXPERIMENTS_name]["merzifon"]["Reads"], 118 )
        self.assertEqual( ribo_info[EXPERIMENTS_name]["ankara"]["Reads"], 12 )
        
        self.assertEqual( ribo_info["Min Read Length"], 2 )
        self.assertEqual( ribo_info["Max Read Length"], 5 )
        self.assertEqual( ribo_info["Reference"], "hg38" )
        
    def test_ribo_attributes(self):
        experiments = self.sample_ribo.experiments
        self.assertEqual(len(experiments), 2)
        
        self.assertTrue("ankara" in experiments)
        self.assertTrue("merzifon" in experiments)
        
        self.assertEqual(self.sample_ribo.reference_name, "hg38")
        
        # Format and version values are changing
        # so we just check we can retrieve them
        # we do NOT check the correctness of the values
        self.assertEqual(type( float(self.sample_ribo.format_version) ), float)
        self.assertEqual(type( self.sample_ribo.ribopy_version), str)
        
        self.assertEqual(self.sample_ribo.minimum_length, 2)
        self.assertEqual(self.sample_ribo.maximum_length, 5)
        
        self.assertEqual(self.sample_ribo.metagene_radius, METAGENE_RADIUS)
        
        self.assertEqual(self.sample_ribo.left_span, LEFT_SPAN)
        self.assertEqual(self.sample_ribo.right_span, RIGHT_SPAN)
        
    def test_transcript_names(self):
        expected_transcripts = ["GAPDH", "VEGFA", "MYC"]
        
        ribo_transcripts = self.sample_ribo.transcript_names
        
        for r in expected_transcripts:
            self.assertTrue( r in ribo_transcripts )
        self.assertEqual(len(ribo_transcripts), 3)
            
        # Also make sure that it works for the cached case
        
        ribo_transcripts_cached = self.sample_ribo.transcript_names
        
        for r in expected_transcripts:
            self.assertTrue( r in ribo_transcripts_cached )
        self.assertEqual(len(ribo_transcripts), 3)
    
    def test_print_info(self):
        import re
        info_str = self.sample_ribo.print_info(return_str = True)
        
        k = re.search(r"Min Read Length[\s]*:[\s]2", info_str)
        self.assertTrue( "ankara" in info_str)
        self.assertTrue( "merzifon" in info_str)
        
        str_buffer = StringIO()
        original_stdout = sys.stdout
        sys.stdout = str_buffer
        self.sample_ribo.print_info()
        captured_output = str_buffer.getvalue()
        sys.stdout = original_stdout
        self.assertTrue("Experiments" in captured_output )
        
    def test_transcript_lengths(self):        
        self.assertEqual( self.sample_ribo.transcript_lengths["GAPDH"], 20 )
        self.assertEqual( self.sample_ribo.transcript_lengths["VEGFA"], 22 )
        self.assertEqual( self.sample_ribo.transcript_lengths["MYC"], 17 )
        
    def test_transcript_index(self):        
        self.assertEqual( self.sample_ribo.transcript_index["GAPDH"], 0 )
        self.assertEqual( self.sample_ribo.transcript_index["VEGFA"], 1 )
        self.assertEqual( self.sample_ribo.transcript_index["MYC"], 2 )    
    
        
        
if __name__ == '__main__':
        
    unittest.main()
