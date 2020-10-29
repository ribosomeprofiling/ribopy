# -*- coding: utf-8 -*-
import unittest

import os
from io import StringIO, BytesIO
import pipes
import shutil

import numpy as np
import pandas as pd
import h5py

from ribopy import create
from ribopy.core.exceptions import NORNASEQ
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.core import create_experiment
from ribopy.settings import *
from ribopy.dump import *
from ribopy.rnaseq import *

from ribopy.merge import merge_ribos

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)

from multilength_test_data import *
import rnaseq_data


###########################################

NPROCESS = 4

class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_len_file     = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file  = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file   = StringIO(READ_SET_1)
        self.alignment_file_2 = StringIO(READ_SET_2)
        self.tsv_1            = StringIO(rnaseq_data.RNASEQ_tsv_1)
        self.tsv_2            = StringIO(rnaseq_data.RNASEQ_tsv_2)
        self.tsv_3            = StringIO(rnaseq_data.RNASEQ_tsv_3)

        self.handle   = h5py.File( BytesIO(), "w" )
        self.handle_2 = h5py.File( BytesIO(), "w" )

        self.rnaseq_reads_handle   = StringIO( rnaseq_data.RNASEQ_READS )
        self.rnaseq_reads_handle_2 = StringIO( rnaseq_data.RNASEQ_READS_2 )

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

        self.ref_len_file.seek(0)
        self.annotation_file.seek(0)

        create.create_ribo(     self.handle_2, "adana",
                                alignment_file  = self.alignment_file_2,
                                reference_name  = "appris_human_v2",
                                length_min      = 2,
                                length_max      = 5,
                                metagene_radius = METAGENE_RADIUS,
                                left_span       = LEFT_SPAN,
                                right_span      = RIGHT_SPAN,
                                lengths_file    = self.ref_len_file,
                                store_coverage  = True,
                                annotation_file = self.annotation_file)

        self.merged_ribo = h5py.File(BytesIO(), "w")
        merge_ribos( self.merged_ribo, [self.handle , self.handle_2] )

    def tearDown(self):
        self.handle.close()

    def initial_rnaseq(self):
        pass


        #set_rnaseq(self.handle, "merzifon", rnaseq_df)

    ### T E S T   R N A S E Q  ################################

    def test_set_rnaseq_from_bed(self):
        set_rnaseq(ribo_handle   = self.handle,
                   name          = "merzifon",
                   rnaseq_reads  = self.rnaseq_reads_handle,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        rnaseq_counts = self.handle[EXPERIMENTS_name]["merzifon"]\
                              [RNASEQ_name][RNASEQ_name][...]
        self.assertTrue( np.isclose( rnaseq_counts[0,:], [1, 0, 2, 1, 3] ).all() )
        self.assertTrue( np.isclose( rnaseq_counts[1,:], [0, 1, 3, 0, 1] ).all() )
        self.assertTrue( np.isclose( rnaseq_counts[2,:], [0, 1, 0, 0, 2] ).all() )

        set_rnaseq(ribo_handle   = self.handle,
                   name          = "merzifon",
                   rnaseq_reads  = self.rnaseq_reads_handle_2,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        rnaseq_counts = self.handle[EXPERIMENTS_name]["merzifon"]\
                              [RNASEQ_name][RNASEQ_name][...]

        self.assertTrue( np.isclose( rnaseq_counts[1,:], [0, 0, 0, 2, 0] ).all() )

    def test_get_rnaseq_from_bed(self):
        set_rnaseq(ribo_handle   = self.handle,
                   name          = "merzifon",
                   rnaseq_reads  = self.rnaseq_reads_handle,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        rnaseq_df = get_rnaseq(self.handle, "merzifon")

        cds_counts = rnaseq_df.loc[ ('merzifon', ), "CDS" ]
        self.assertTrue( np.isclose( cds_counts, [2, 3, 0] ).all() )


    def test_get_rnaseq_from_multiple(self):
        set_rnaseq(ribo_handle   = self.merged_ribo,
                   name          = "merzifon",
                   rnaseq_reads  = self.rnaseq_reads_handle,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        set_rnaseq(ribo_handle   = self.merged_ribo,
                   name          = "adana",
                   rnaseq_reads  = self.rnaseq_reads_handle_2,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        rnaseq_counts =  get_rnaseq(self.merged_ribo)

        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "GAPDH"), "UTR5"] , 1 ) )
        self.assertTrue( np.isclose( rnaseq_counts.loc[("adana", "VEGFA"), "CDS"] , 0 ) )

    def test_delete_rnaseq(self):
        set_rnaseq(ribo_handle   = self.merged_ribo,
                   name          = "merzifon",
                   rnaseq_reads  = self.rnaseq_reads_handle,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        set_rnaseq(ribo_handle   = self.merged_ribo,
                   name          = "adana",
                   rnaseq_reads  = self.rnaseq_reads_handle_2,
                   format        = "bed",
                   rnaseq_counts = None,
                   sep           = "\t")

        rnaseq_counts =  get_rnaseq(self.merged_ribo, "merzifon")
        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "GAPDH"), "UTR3"] , 3 ) )
        delete_rnaseq(self.merged_ribo, "merzifon")

        with self.assertRaises( NORNASEQ ) as context:
            get_rnaseq(self.merged_ribo, "merzifon")

        rnaseq_counts =  get_rnaseq(self.merged_ribo, "adana")
        self.assertTrue( np.isclose(rnaseq_counts.loc[("adana", "MYC"), "CDS"] , 0)  )

    def test_set_rnaseq_from_tsv_1(self):
        set_rnaseq(ribo_handle   = self.handle,
                   rnaseq_reads  = None,
                   name          = "merzifon",
                   rnaseq_counts = self.tsv_1,
                   format        = "bed",
                   sep           = "\t")

        rnaseq_counts =  get_rnaseq(self.handle, "merzifon")
        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "GAPDH"), "UTR5_junction"] , 7.2 ) )
        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "MYC"), "UTR5"] , 0 ) )


    def test_set_rnaseq_from_tsv_2(self):
        set_rnaseq(ribo_handle   = self.handle,
                   rnaseq_reads  = None,
                   name          = "merzifon",
                   rnaseq_counts = self.tsv_2,
                   format        = "bed",
                   sep           = "\t")

        rnaseq_counts =  get_rnaseq(self.handle, "merzifon")
        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "VEGFA"), "CDS"] , 5.7 ) )
        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "MYC"), "CDS"] , 0 ) )
        self.assertTrue( np.isclose( rnaseq_counts.loc[("merzifon", "GAPDH"), "UTR5"] , 10 ) )


    def test_set_rnaseq_from_tsv_3(self):
        """
        Make sure we reject tables with missing columns.
        """
        with self.assertRaises( RiboBaseError ) as context:
            set_rnaseq(ribo_handle   = self.handle,
                       rnaseq_reads  = None,
                       name          = "merzifon",
                       rnaseq_counts = self.tsv_3,
                       format        = "bed",
                       sep           = "\t")

    def test_rnaseq_from_bam(self):
        transcript_lengths_file = "trans_len.tsv"
        bam_file                = "rnaseq.bam"

        with open(transcript_lengths_file , "w") as output_stream:
            print(TRANSCRIPT_LENGTHS,
                  file = output_stream, end = "")

        if not shutil.which('bamToBed'):
            raise FileNotFoundError("Could not find the executable bamToBed")
        if not shutil.which('bedToBam'):
            raise FileNotFoundError("Could not find the executable bedToBam")

        pipe_t = pipes.Template()
        pipe_t.append("bedToBam -g {}"\
                      .format( transcript_lengths_file), "--")

        f = pipe_t.open(bam_file, 'w')
        f.write(rnaseq_data.RNASEQ_READS)
        f.close()

        set_rnaseq(ribo_handle   = self.handle,
                   rnaseq_reads  = bam_file,
                   name          = "merzifon",
                   rnaseq_counts = None,
                   format        = "bam")

        rnaseq_counts =  self.handle[EXPERIMENTS_name]["merzifon"]\
                              [RNASEQ_name][RNASEQ_name][...]

        self.assertTrue( np.isclose( rnaseq_counts[0,:], [1, 0, 2, 1, 3] ).all() )
        self.assertTrue( np.isclose( rnaseq_counts[1,:], [0, 1, 3, 0, 1] ).all() )
        self.assertTrue( np.isclose( rnaseq_counts[2,:], [0, 1, 0, 0, 2] ).all() )

        for f in (transcript_lengths_file, bam_file):
            os.remove(f)




if __name__ == '__main__':

    unittest.main()
