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
from ribopy.core.get_gadgets import *
from ribopy.metadata import set_metadata
from ribopy.create import set_annotation, \
                          check_annotation, \
                          post_check_annotation, \
                          set_reference_names_and_lengths, \
                          initialize

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from annotation_test_data import *
from api_test_base import ApiTestBase

"""
Now we allow empty 3' UTRs in our annotation.
This required modifying our annotration validation and quantification steps a little.
So we added additional test where our system can correctly handle genes with
"""

####################################################################
"""
class ApiTestBase(unittest.TestCase):

    def setUp(self):
        self.tmp_files = list()


        self.ref_len_file       = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file    = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file_1   = StringIO(READ_SET_1)
        self.alignment_file_2   = StringIO(READ_SET_2)

        self.handle   = h5py.File(BytesIO(), "w" )
        self.handle_2 = h5py.File(BytesIO(), "w" )

        create.create_ribo(
                ribo            = self.handle,
                experiment_name = "Experiment_1",
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

        create.create_ribo(
                ribo            = self.handle_2,
                experiment_name = "Experiment_2",
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


        self.merged_ribo.close()

        self.sample_ribo = Ribo( self.merged_io )

    def tearDown(self):
        self.handle.close()
        self.handle_2.close()

    def test_check_annotation(self):
        print("Generating annotation...")

        h5_handle        = h5py.File(BytesIO(), "w" )

        initialize(h5_handle, experiment_name = "dummy",
                              reference_name  = "test_ref")

        (ref_names, ref_lengths) = \
           set_reference_names_and_lengths(h5_handle , self.ref_len_file)
        self.ref_len_file.seek(0)

        annotation_lines = TRANSCRIPT_ANNOTATION.split("\n")
        # This should pass
        check_annotation( h5_handle, annotation_lines )

        #with self.assertRaises(AnnotationError):

"""

"""
class BasicAnnotationTest(unittest.TestCase):
    def setUp(self):

        self.tmp_files = list()


        self.ref_len_file       = StringIO(SINGLE_TRANSCRIPT_LENGTHS)
        self.annotation_file    = StringIO(SINGLE_TRANSCRIPT_ANNOTATION)

        self.h5_handle  = h5py.File(BytesIO(), "w" )

        initialize(self.h5_handle, experiment_name = "dummy",
                              reference_name  = "test_ref")

        (ref_names, ref_lengths) = \
           set_reference_names_and_lengths(self.h5_handle , self.ref_len_file)
        self.ref_len_file.seek(0)

        annotation_lines = SINGLE_TRANSCRIPT_ANNOTATION.split("\n")

    def test_check_annotation(self):

        annotation_lines = SINGLE_TRANSCRIPT_ANNOTATION.split("\n")
        # This should pass
        check_annotation( self.h5_handle, annotation_lines )

        annotation_lines_2 = SINGLE_TRANSCRIPT_ANNOTATION_2.split("\n")
        # This should pass
        with self.assertRaises(AnnotationError):
            check_annotation( self.h5_handle, annotation_lines_2 )

        annotation_lines = NOUTR5_SINGLE_TRANSCRIPT_ANNOTATION.split("\n")
        # This should pass
        check_annotation( self.h5_handle, annotation_lines )

        annotation_lines = NOUTR3_SINGLE_TRANSCRIPT_ANNOTATION.split("\n")
        # This should pass
        check_annotation( self.h5_handle, annotation_lines )
"""

class GenericAnnotationTest(unittest.TestCase):
    def setUp(self):

        self.tmp_files = list()


        self.ref_len_file       = StringIO(GENERIC_TRANSCRIPT_LENGTHS)
        self.annotation_file    = StringIO(GENERIC_TRANSCRIPT_ANNOTATION)
        self.alignment_file_1   = StringIO(READ_SET_1)
        self.alignment_file_2   = StringIO(READ_SET_2)

        self.handle_io  = BytesIO()
        self.handle     = h5py.File(self.handle_io, "w" )
        self.handle_2   = h5py.File(BytesIO(),      "w" )
        self.h5_handle  = h5py.File(BytesIO(),      "w" )

        initialize(self.h5_handle, experiment_name = "dummy",
                              reference_name  = "test_ref")

        (ref_names, ref_lengths) = \
           set_reference_names_and_lengths(self.h5_handle , self.ref_len_file)
        self.ref_len_file.seek(0)

        self.annotation_lines = GENERIC_TRANSCRIPT_ANNOTATION.split("\n")

        create.create_ribo(
                ribo            = self.handle,
                experiment_name = "experiment-1",
                alignment_file  = self.alignment_file_1,
                reference_name  = "hg38",
                lengths_file    = self.ref_len_file,
                annotation_file = self.annotation_file,
                metagene_radius = METAGENE_RADIUS,
                left_span       = LEFT_SPAN,
                right_span      = RIGHT_SPAN,
                length_min      = LENGTH_MIN,
                length_max      = LENGTH_MAX,
                store_coverage  = True,
                nprocess        = NPROCESS,
                tmp_file_prefix = "")

        self.handle.close()
        self.sample_ribo = Ribo(self.handle_io)

        # Get the  region counts

        self.cds_counts = \
           self.sample_ribo.get_region_counts( region_name    = CDS_name,
                                               sum_references = False)

        self.utr5_counts = \
           self.sample_ribo.get_region_counts( region_name    = UTR5_name,
                                               sum_references = False)

        self.utr3_counts = \
           self.sample_ribo.get_region_counts( region_name    = UTR3_name,
                                               sum_references = False)

        self.utr5j_counts = \
           self.sample_ribo.get_region_counts( region_name    = UTR5_JUNCTION_name,
                                               sum_references = False)

        self.utr3j_counts = \
           self.sample_ribo.get_region_counts( region_name    = UTR3_JUNCTION_name,
                                               sum_references = False)

    def test_check_annotation(self):

        # This should pass
        check_annotation( self.h5_handle, self.annotation_lines)


    def test_check_no_CDS(self):

        self.assertEqual(self.cds_counts.loc["Gene_Only-CDS"][0], 5)

        self.assertEqual(self.utr5_counts.loc["Gene_Only-CDS"][0], 0)

        self.assertEqual(self.utr5j_counts.loc["Gene_Only-CDS"][0], 4)

        self.assertEqual(self.utr3j_counts.loc["Gene_Only-CDS"][0], 1)

        self.assertEqual(self.utr3_counts.loc["Gene_Only-CDS"][0], 0)


    def test_NO_UTR5(self):
        self.assertEqual(self.cds_counts.loc["Gene_NO-UTR5"][0], 1)

        self.assertEqual(self.utr5_counts.loc["Gene_NO-UTR5"][0], 0)

        self.assertEqual(self.utr5j_counts.loc["Gene_NO-UTR5"][0], 2)

        self.assertEqual(self.utr3j_counts.loc["Gene_NO-UTR5"][0], 0)

        self.assertEqual(self.utr3_counts.loc["Gene_NO-UTR5"][0], 1)


    def test_NO_UTR3(self):
        self.assertEqual(self.cds_counts.loc["Gene_NO-UTR3"][0], 2)

        self.assertEqual(self.utr5_counts.loc["Gene_NO-UTR3"][0], 2)

        self.assertEqual(self.utr5j_counts.loc["Gene_NO-UTR3"][0], 3)

        self.assertEqual(self.utr3j_counts.loc["Gene_NO-UTR3"][0], 1)

        self.assertEqual(self.utr3_counts.loc["Gene_NO-UTR3"][0], 0)

    def test_short_UTR3(self):
        self.assertEqual(self.cds_counts.loc["Gene_short_UTR3"][0], 0)

        self.assertEqual(self.utr5_counts.loc["Gene_short_UTR3"][0], 0)

        self.assertEqual(self.utr5j_counts.loc["Gene_short_UTR3"][0], 0)

        self.assertEqual(self.utr3j_counts.loc["Gene_short_UTR3"][0], 1)

        self.assertEqual(self.utr3_counts.loc["Gene_short_UTR3"][0], 0)

if __name__ == '__main__':
    unittest.main()
