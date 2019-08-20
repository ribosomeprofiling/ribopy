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
from ribopy.core import create_experiment
from ribopy.settings import *
from ribopy.metadata import *

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

        self.handle = h5py.File(BytesIO())

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

    def test_set_metadata(self):
        exp_handle  = self.handle[EXPERIMENTS_name]["merzifon"]
        ribo_handle = self.handle

        metadata_1  = \
            "cell-count : 5\nscientist : john doe\nsequencer : MiSeq"
        set_metadata(self.handle, "merzifon", metadata_1)

        self.assertTrue( "MiSeq" in exp_handle.attrs[USER_METADATA])

        metadata_2_str = \
            "cell-count : 10\nscientist : who doe\nsequencer : HiSeq-4000"
        metadata_2 = StringIO(metadata_2_str)
        set_metadata(self.handle, "merzifon", metadata_2)

        self.assertTrue( "MiSeq" not in exp_handle.attrs[USER_METADATA])

        self.assertTrue( "HiSeq-4000" in exp_handle.attrs[USER_METADATA])

        ribo_metadata = "pipelineversion : v1.0.2\ndate: 3-5-2018"

        set_metadata(self.handle, "", ribo_metadata)
        self.assertTrue("pipelineversion" in ribo_handle.attrs["metadata"])

        # Invalid Metadata

        invalid_metadata = "nonsens : 3\tsome-sense : 89"
        error = set_metadata(self.handle, "", invalid_metadata) 
        self.assertTrue("Invalid" in error)

    def test_get_metadata(self):
        exp_handle  = self.handle[EXPERIMENTS_name]["merzifon"]
        ribo_handle = self.handle

        metadata_1  = \
            "cell-count : 5\nscientist : john doe\nsequencer : MiSeq"
        set_metadata(self.handle, "merzifon", metadata_1)

        metadata_back_1 = get_metadata( self.handle, "merzifon" )
        self.assertEqual( str(metadata_back_1["cell-count"]), "5")
        self.assertEqual(metadata_back_1["sequencer"], "MiSeq")

        ribo_metadata   = "pipelineversion : v1.0.2\ndate: 3-5-2018"

        set_metadata(self.handle, "", ribo_metadata)
        ribo_metadata_back = get_metadata(self.handle, "")
        self.assertEqual(ribo_metadata_back["date"], "3-5-2018")

    def test_get_empty_metadata(self):
        exp_handle  = self.handle[EXPERIMENTS_name]["merzifon"]
        ribo_handle = self.handle

        metadata_back_1 = get_metadata( self.handle, "merzifon" )

        self.assertEqual(metadata_back_1, None)




    

if __name__ == '__main__':
        
    unittest.main()
