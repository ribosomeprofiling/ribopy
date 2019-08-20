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
from ribopy.metadata import set_metadata

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from multilength_test_data import *

NPROCESS = 1

class ApiTestBase(unittest.TestCase):

    def setUp(self):
        self.tmp_files = list()


        self.ref_len_file       = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file    = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file_1   = StringIO(READ_SET_1)
        self.alignment_file_2   = StringIO(READ_SET_2)

        self.handle   = h5py.File(BytesIO() )
        self.handle_2 = h5py.File(BytesIO() )

        create.create_ribo(
                ribo            = self.handle, 
                experiment_name = "merzifon", 
                alignment_file  = self.alignment_file_1,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file, 
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS, 
                left_span       = LEFT_SPAN, 
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                store_coverage  = False,  
                nprocess        = NPROCESS,
                tmp_file_prefix = "")
        self.ref_len_file.seek(0)
        self.annotation_file.seek(0)
        self.alignment_file_1 = StringIO(READ_SET_1)
        
        set_metadata(self.handle, "merzifon", METADATA_EXPERIMENT_DICT)
        
        create.create_ribo(
                ribo            = self.handle_2, 
                experiment_name = "ankara", 
                alignment_file  = self.alignment_file_2,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file, 
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS, 
                left_span       = LEFT_SPAN, 
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                store_coverage  = True,
                nprocess        = NPROCESS,
                tmp_file_prefix = "")
        self.ref_len_file.seek(0)
        self.annotation_file.seek(0)
        self.alignment_file_2 = StringIO(READ_SET_2)

        self.merged_io   = BytesIO()
        self.merged_ribo = h5py.File(self.merged_io, "w")
        merge_ribos( self.merged_ribo, [self.handle , self.handle_2] )
        set_metadata(self.merged_ribo, 
                     name     = None, 
                     metadata = RIBO_METADATA_STR_1)
        self.merged_ribo.close()
        
        self.sample_ribo = Ribo( self.merged_io )

    def tearDown(self):
        self.handle.close()
        self.handle_2.close()
