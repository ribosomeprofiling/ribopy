# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO, BytesIO

import numpy as np
import h5py

from ribopy import create
from ribopy.merge import merge_ribos
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries,\
                                         get_experiment_names
from ribopy.core import create_experiment
from ribopy.settings import *
from ribopy.dump import *

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)

from multilength_test_data import *


###########################################

INCOMPATIBLE_TRANSCRIPT_LENGTHS=\
"""GAPDH    20
VEGFA   22
MYC 17
EXTRA   56"""

incomp_len_file = "incompatible_lengths.tsv"

INCOMPATIBLE_TRANSCRIPT_ANNOTATION=\
"""GAPDH    0   5  UTR5    0   +
GAPDH   5  15    CDS 0   +
GAPDH   15  20    UTR3 0   +
VEGFA    0   4  UTR5    0   +
VEGFA   4  16    CDS 0   +
VEGFA   16  22    UTR3 0   +
MYC    0   3  UTR5    0   +
MYC   3  10    CDS 0   +
MYC   10  17    UTR3 0   +
EXTRA   0  20    UTR5 0   +
EXTRA   20  49    CDS 0   +
EXTRA   49  56    UTR3 0   +"""

incompatible_annotation = "incomp_annotation.bed"

###########################################

NPROCESS = 4

class TestCreate(unittest.TestCase):

    def setUp(self):
        self.tmp_files = list()


        self.ref_len_file       = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file    = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file_1   = StringIO(READ_SET_1)
        self.alignment_file_2   = StringIO(READ_SET_2)
        self.in_ref_len_file    = StringIO(INCOMPATIBLE_TRANSCRIPT_LENGTHS)
        self.in_annotation_file = StringIO(INCOMPATIBLE_TRANSCRIPT_ANNOTATION)

        self.handle   = h5py.File(BytesIO(), "w" )
        self.handle_2 = h5py.File(BytesIO(), "w" )

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
                nprocess        = NPROCESS,
                tmp_file_prefix = "")
        self.ref_len_file.seek(0)
        self.annotation_file.seek(0)
        self.alignment_file_1 = StringIO(READ_SET_1)

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
                nprocess        = NPROCESS,
                tmp_file_prefix = "")
        self.ref_len_file.seek(0)
        self.annotation_file.seek(0)
        self.alignment_file_2 = StringIO(READ_SET_2)

        self.merged_ribo = h5py.File(BytesIO(), "w")
        merge_ribos( self.merged_ribo, [self.handle , self.handle_2] )


    def tearDown(self):
        self.handle.close()
        self.handle_2.close()
        self.merged_ribo.close()

    def test_merge_libnames(self):
        """
        Makes sure that the merged ribo file has the
        expected experiment names
        """

        self.handle.close()
        self.handle_2.close()

        merged_exp_names = set(get_experiment_names(self.merged_ribo)  )
        expected_libnames = set( ("ankara", "merzifon") )

        self.assertTrue( merged_exp_names == expected_libnames )


    def test_merge_regioncounts(self):
        cds_counts = dump_region(
                        ribo_handle     = self.merged_ribo,
                        region_name     = CDS_name,
                        sum_lengths     = True, sum_references = False,
                        range_lower     = 2, range_upper =5,
                        experiment_list = [])

        self.assertListEqual( list(cds_counts["ankara"]) , [5 , 3, 0])
        self.assertListEqual( list(cds_counts["merzifon"]) , [13 , 10, 3])

    def test_merge_start_site(self):

        start_site_cov_dict = dump_metagene(
                                   ribo_handle     = self.merged_ribo,
                                   site_type       = REF_DG_START_SITE_COV,
                                   sum_lengths     = True,
                                   sum_references  = True,
                                   range_lower     = 3,
                                   range_upper     = 3,
                                   experiment_list = [] )

        self.assertListEqual(
             list(start_site_cov_dict["ankara"].iloc[0]), [0,0,0,0,0] )

        self.assertListEqual(
            list(start_site_cov_dict["merzifon"].iloc[0]), [3,5,3,2,0])

    def test_merge_single_ribo(self):
        single_merged = h5py.File (BytesIO() , "w" )

        merge_ribos( single_merged, [self.handle] )

        self.assertTrue( "merzifon" in get_experiment_names(single_merged) )


    def test_merge_incompatible_length(self):
        handle_incompatible = h5py.File(BytesIO(),"w" )
        ribo_out            = h5py.File(BytesIO(), "w" )

        create.create_ribo(
                ribo            = handle_incompatible,
                experiment_name = "izmir",
                alignment_file  = self.alignment_file_2,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file,
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN,
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 3, # incompatible parameter
                nprocess        = NPROCESS,
                tmp_file_prefix = "")

        #merge_ribos( ribo_out, [self.handle, self.handle_2, handle_incompatible] )

        with self.assertRaises(ValueError):
            merge_ribos(ribo_out,
                        [self.handle, self.handle_2, handle_incompatible] )


    def test_merge_incompatible_samelib(self):
        handle_incompatible = h5py.File(BytesIO(), "w" )
        ribo_out            = h5py.File(BytesIO(), "w" )

        create.create_ribo(
                ribo            = handle_incompatible,
                experiment_name    = "ankara", # exception reason
                alignment_file  = self.alignment_file_2,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file,
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN,
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                nprocess        = NPROCESS,
                tmp_file_prefix = "")

        #merge_ribos( ribo_out, [self.handle, self.handle_2, handle_incompatible] )

        with self.assertRaises(ValueError):
            merge_ribos( ribo_out, [self.handle, self.handle_2,
                                    handle_incompatible] )



    def test_merge_incompatible_refname(self):
        handle_incompatible = h5py.File(BytesIO(), "w" )
        ribo_out            = h5py.File(BytesIO(), "w" )

        create.create_ribo(
                ribo            = handle_incompatible,
                experiment_name = "izmir",
                alignment_file  = self.alignment_file_2,
                reference_name  = "hg19", # incompatible reference name
                lengths_file    = self.ref_len_file,
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN,
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                nprocess        = NPROCESS,
                tmp_file_prefix = "")


        with self.assertRaises(ValueError):
            merge_ribos(ribo_out, [self.handle,
                                   self.handle_2, handle_incompatible] )


    def test_merge_incompatible_leftspan(self):
        handle_incompatible = h5py.File(BytesIO(), "w" )
        ribo_out            = h5py.File(BytesIO(), "w" )

        create.create_ribo(
                ribo            = handle_incompatible,
                experiment_name = "izmir",
                alignment_file  = self.alignment_file_2,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file,
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN + 3, # incomp span
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                nprocess        = NPROCESS,
                tmp_file_prefix = "")


        with self.assertRaises(ValueError):
            merge_ribos(ribo_out, [self.handle,
                                   self.handle_2, handle_incompatible] )



    def test_merge_incompatible_reference(self):
        handle_incompatible = h5py.File(BytesIO(), "w" )
        ribo_out            = h5py.File(BytesIO(), "w" )

        create.create_ribo(
                ribo            = handle_incompatible,
                experiment_name = "izmir",
                alignment_file  = self.alignment_file_2,
                reference_name  = "hg38",
                lengths_file    = self.in_ref_len_file, # incomp lens
                annotation_file = self.in_annotation_file, # incomp annot.
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN,
                right_span      = RIGHT_SPAN,
                length_min      = 2,
                length_max      = 5,
                nprocess        = NPROCESS,
                tmp_file_prefix = "")


        with self.assertRaises(ValueError):
            merge_ribos(ribo_out, [self.handle,
                                   self.handle_2, handle_incompatible] )



if __name__ == '__main__':

    unittest.main()
