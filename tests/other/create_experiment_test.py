# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO, BytesIO

import yaml
import numpy as np
import h5py

from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *
from ribopy.rnaseq import get_rnaseq

import sys
test_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir)

from multilength_test_data import *
from ribopy.core import create_experiment

###########################################

NPROCESS = 4

class TestCreate(unittest.TestCase):

    def setUp(self):
        self.ref_length_file  = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file  = StringIO(TRANSCRIPT_ANNOTATION)
        self.annotation_lines = TRANSCRIPT_ANNOTATION.split("\n")
        self.alignment_file   = StringIO(READ_SET_1)
        self.metadata_file    = StringIO(METADATA_EXPERIMENT_STR_1)

        self.handle = h5py.File(BytesIO(), "w")
        create.initialize( self.handle, "merzifon", "hg38" )
        create.set_reference_names_and_lengths(
                 self.handle , self.ref_length_file)
        create.set_annotation( self.handle, self.annotation_lines )

        ref_names          = get_reference_names(self.handle)
        ref_lengths        = get_reference_lengths(self.handle)
        region_coordinates = get_region_boundaries(self.handle)

        exp_handle = self.handle[EXPERIMENTS_name]["merzifon"]

        create_experiment.create_experiment(
            ribo_exp_handle       = exp_handle,
            experiment_name       = "merzifon",
            alignment_file_handle = self.alignment_file,
            ref_names             = ref_names,
            ref_lengths           = ref_lengths,
            region_coordinates    = region_coordinates,
            metagene_radius       = METAGENE_RADIUS,
            left_span             = LEFT_SPAN,
            right_span            = RIGHT_SPAN,
            length_min            = 2,
            length_max            = 5,
            metadata              = self.metadata_file,
            store_coverage        = True,
            nprocess              = 4)

    def tearDown(self):
        self.handle.close()

    def test_start_sites(self):
        start_site_0       = np.array([0,1,0,1,2,0,4,3,0,1,0,0])
        start_site_minus_2 = np.array([1,0,0,3,0,0,2,0,0,1,0,0])
        start_site_plus_2  = np.array([1,0,0,0,0,0,1,0,0,0,0,0])

        metagene_handle = self.handle[EXPERIMENTS_name]["merzifon"][METAGENE_name]
        start_sites     = metagene_handle[REF_DG_START_SITE_COV]

        self.assertTrue( np.all(start_site_0       == start_sites[:,2] ) )
        self.assertTrue( np.all(start_site_minus_2 == start_sites[:,0] ) )
        self.assertTrue( np.all(start_site_plus_2  == start_sites[:,4] ) )

    def test_stop_site(self):
        stop_site_0        = np.array([1,1,2,2,2,1,2,0,3,0,2,7])
        stop_site_minus_2  = np.array([2,0,0,4,0,0,4,0,0,0,0,0])
        stop_site_plus_2   = np.array([0,0,0,0,0,0,0,0,0,0,0,0])

        metagene_handle = self.handle[EXPERIMENTS_name]["merzifon"][METAGENE_name]
        stop_sites      = metagene_handle[REF_DG_STOP_SITE_COV]

        self.assertTrue( np.all(stop_site_0       == stop_sites[:,2] ) )
        self.assertTrue( np.all(stop_site_minus_2 == stop_sites[:,0] ) )
        self.assertTrue( np.all(stop_site_plus_2  == stop_sites[:,4] ) )

    def test_region_counts(self):
        region_counts_handle = \
          self.handle[EXPERIMENTS_name]["merzifon"][REF_DG_REGION_COUNTS]
        region_counts_all = region_counts_handle[REF_DG_REGION_COUNTS][...]

        len_3_offset = 3
        #( 2, 7, 3, 9 , 0 )
        self.assertTrue(np.all(GAPDH_counts_length_3 == \
                         region_counts_all[len_3_offset, :] ) )

        #( 0, 6, 2, 2, 0 )
        self.assertTrue( np.all( VEGFA_counts_length_3 == \
                          region_counts_all[len_3_offset +1, :] ))

        # ( 0, 1, 1, 1, 0 )
        self.assertTrue(np.all(MYC_counts_length_3 ==\
                          region_counts_all[len_3_offset +2, :]) )


    def test_coverage(self):
        """
        TRANSCRIPT_LENGTHS=
        GAPDH    20
        VEGFA   22
        MYC 17

        Length_range = [2, 5]
        """
        coverage_handle = \
            self.handle[EXPERIMENTS_name]["merzifon"][REF_DG_COVERAGE]

        all_coverage = coverage_handle[REF_DG_COVERAGE][...]

        transcriptome_size = 20 + 22 + 17

        # GAPDH, len = 2, pos = 4
        # offset is 0 brcause 2 is the first length
        # and gapdh is the first gene
        offset   = 0
        coverage = all_coverage[offset + 4]
        self.assertEqual(coverage, 2)

        # GAPDH, len = 4, pos = 7
        offset   = transcriptome_size * 2
        coverage = all_coverage[offset + 7]
        self.assertEqual(coverage, 1)

        # VEGF, len = 5, pos = 10
        offset   = transcriptome_size * 3 + 20
        coverage = all_coverage[offset + 10]
        self.assertEqual(coverage, 4)

        # VEGF, len = 5, pos = 9
        offset   = transcriptome_size * 3 + 20
        coverage = all_coverage[offset + 9]
        self.assertEqual(coverage, 0)

        # VEGF, len = 5, pos = 11
        offset   = transcriptome_size * 3 + 20
        coverage = all_coverage[offset + 11]
        self.assertEqual(coverage, 0)

    def test_total_reads(self):
        exp_handle  = self.handle[EXPERIMENTS_name]["merzifon"]
        total_reads = exp_handle.attrs[ATTRS_TOTAL_READS]

        # In the test data there are 128 lines
        self.assertEqual(total_reads, 118)


    def test_exp_metadata(self):
        exp_handle  = self.handle[EXPERIMENTS_name]["merzifon"]

        metadata = yaml.safe_load(exp_handle.attrs.get(USER_METADATA, ""))

        self.assertEqual( metadata["cell_line"], "HeLa" )
        self.assertEqual( metadata["link"],
                            "https://www.encodeproject.org/" )
        self.assertEqual( metadata["digestion_duration"], "5 min" )

        total_reads = exp_handle.attrs[ATTRS_TOTAL_READS]
        self.assertEqual(total_reads, 118)


###################################################################################


class TestCreateFromMetaFile(unittest.TestCase):

    def setUp(self):
        self.annotation_lines = TRANSCRIPT_ANNOTATION.split("\n")
        self.ref_length_file  = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file  = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file   = StringIO(READ_SET_1)
        self.metadata_file    = StringIO(METADATA_EXPERIMENT_STR_1)

        self.handle = h5py.File(BytesIO(), "w")
        create.initialize( self.handle, "merzifon", "hg38" )
        create.set_reference_names_and_lengths(self.handle ,
                                               self.ref_length_file)
        create.set_annotation( self.handle, self.annotation_lines )

        ref_names          = get_reference_names(self.handle)
        ref_lengths        = get_reference_lengths(self.handle)
        region_coordinates = get_region_boundaries(self.handle)

        exp_handle = self.handle[EXPERIMENTS_name]["merzifon"]

        create_experiment.create_experiment(
            ribo_exp_handle       = exp_handle,
            experiment_name          = "merzifon",
            alignment_file_handle = self.alignment_file,
            ref_names             = ref_names,
            ref_lengths           = ref_lengths,
            region_coordinates    = region_coordinates,
            metagene_radius       = METAGENE_RADIUS,
            left_span             = LEFT_SPAN,
            right_span            = RIGHT_SPAN,
            length_min            = 2,
            length_max            = 5,
            metadata              = self.metadata_file,
            nprocess              = 4)


    def tearDown(self):
        self.handle.close()

    def test_exp_metadata(self):
        metadata_str = self.handle[EXPERIMENTS_name]["merzifon"].attrs[USER_METADATA]
        metadata     = yaml.safe_load(metadata_str)

        self.assertEqual( metadata["cell_line"], "HeLa" )
        self.assertEqual( metadata["digestion_enzyme"], "HindIII" )
        self.assertEqual( metadata["digestion_duration"], "5 min" )

        """
        self.assertEqual( metadata_handle["cell_line"] , "HeLa" )
        self.assertEqual( metadata_handle["digestion_enzyme"], "HindIII" )
        self.assertEqual( metadata_handle["digestion_duration"], "5 min" )
        self.assertEqual( metadata_handle["link"],
                          "https://www.encodeproject.org/" )
        """




if __name__ == '__main__':

    unittest.main()
