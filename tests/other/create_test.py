# -*- coding: utf-8 -*-
import unittest
from io import StringIO, BytesIO
import os


import pandas as pd
import numpy as np
import h5py


from ribopy import create
from ribopy.settings import *
from ribopy.core.exceptions import *

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)

from test_data import *
###########################################################

len_min = 5
len_max = 35

class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_len_file     = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file  = StringIO(TRANSCRIPT_ANNOTATION)
        self.annotation_lines = StringIO(TRANSCRIPT_ANNOTATION).readlines()
        self.alignment_file   = StringIO(ALIGNMENT_FILE_1)
        self.handle           = h5py.File(BytesIO())
        
    def tearDown(self):
        self.handle.close()


    def test_init(self):
        create.initialize(self.handle, "merzifon", 
                          reference_name = "appris_human_v2")
        self.assertIn("merzifon", 
                      self.handle[EXPERIMENTS_name].keys())
        self.assertIn(REFERENCE_name, self.handle.keys())
        
    def test_set_reference_names_and_lengths(self):
        create.initialize(self.handle, "merzifon", 
                          reference_name = "appris_human_v2")
        create.set_reference_names_and_lengths(self.handle, self.ref_len_file)

        ref_h     = self.handle[REFERENCE_name]
        ref_names = ref_h[REF_DG_REFERENCE_NAMES][...].astype(str) 

        self.assertIn( "GAPDH",  ref_names[0])
        self.assertIn( "BRCA",  ref_names[3])
        self.assertEqual(4, len(ref_names) )

        ref_lengths = ref_h[REF_DG_REFERENCE_LENGTHS][...]
        self.assertEqual(875,     ref_lengths[1])
        self.assertEqual(1462,    ref_lengths[2])
        self.assertNotEqual(1462, ref_lengths[0])

    def test_create_ribo_file(self):
        test_ribo_file = "test_create_file.ribo"
        create.create_ribo_file(test_ribo_file, "iconia", 
                                alignment_file  = self.alignment_file,
                                reference_name  = "appris_human_v2",
                                length_min      = len_min,
                                length_max      = len_max,
                                metagene_radius = 3,
                                left_span       = 3, 
                                right_span      = 2,
                                lengths_file    = self.ref_len_file,
                                annotation_file = self.annotation_file)
        try:
            with h5py.File(test_ribo_file) as h:
                self.assertEqual( tuple(h[EXPERIMENTS_name].keys())[0] , "iconia")
        finally:
            os.remove(test_ribo_file)


    def test_get_reference_names(self):
        create.initialize(self.handle, "merzifon", 
                          reference_name = "appris_human_v2")
        create.set_reference_names_and_lengths(self.handle, self.ref_len_file)

        create.set_annotation( h5_handle        = self.handle, 
                               annotation_lines = self.annotation_lines )

        ref_array = self.handle[REFERENCE_name][REF_ANNOTATION_NAME][...]
        self.assertEqual( ref_array[0,0], 50 )
        self.assertEqual( ref_array[1,0], 0 )
        self.assertEqual( ref_array[1,2], 875 )
        self.assertEqual( ref_array[2,1], 1251 )
        self.assertEqual( ref_array[3,2], 565 )


    def test_set_coverage_vectors(self):

        create.initialize(self.handle, "merzifon", 
                          reference_name = "appris_human_v2")

        create.set_reference_names_and_lengths(self.handle, self.ref_len_file)
        create.set_annotation( h5_handle        = self.handle, 
                               annotation_lines = self.annotation_lines )

        create.set_coverage_vectors(self.handle, 5)

        start_site_cov = self.handle[REFERENCE_name][REF_DG_START_SITE_COV][...]
        stop_site_cov  = self.handle[REFERENCE_name][REF_DG_STOP_SITE_COV][...]
        self.assertListEqual( list(start_site_cov) , [2,2,2,3,3,4,4,4,4,4,4] )
        self.assertListEqual( list(stop_site_cov), [4,4,4,4,4,4,4,4,4,4,3] )
        self.assertTrue( list(stop_site_cov) != [6,4,4,4,4,4,4,4,4,4,3] )
        

    def test_invalid_name_check(self):
        test_ribo_file = "test_create_file.ribo"

        with self.assertRaises(RiboBaseError) as exception_context:

            create.create_ribo_file(test_ribo_file, 
                                    experiment_name = "this,that", 
                                    alignment_file  = self.alignment_file,
                                    reference_name  = "appris_human_v2",
                                    length_min      = len_min,
                                    length_max      = len_max,
                                    metagene_radius = 3,
                                    left_span       = 3, 
                                    right_span      = 2,
                                    lengths_file    = self.ref_len_file,
                                    annotation_file = self.annotation_file)

        if os.path.isfile(test_ribo_file):
            os.remove(test_ribo_file)

if __name__ == '__main__':
        
    unittest.main()
