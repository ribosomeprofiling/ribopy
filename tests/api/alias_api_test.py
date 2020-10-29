# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO, BytesIO

import h5py
import pandas as pd

from ribopy import Ribo
from ribopy import create
from ribopy.merge import merge_ribos
from ribopy.settings import *
from ribopy.api.alias import ReferenceAlias
from ribopy.core.exceptions import *
from ribopy.metadata import set_metadata

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)

from multilength_test_data import *
from api_test_base import ApiTestBase

####################################################################

NPROCESS = 1


class TestAliasClass(unittest.TestCase):
    """
    We only the Alias class in this Test Class
    """
    def setUp(self):
        self.t_names = tuple(map(lambda x: x.split("\t")[0] ,
                                  TRANSCRIPT_LENGTHS.split("\n") ))

    def test_lower_case(self):

        t_lower = lambda x: x.lower()

        ref_alias = ReferenceAlias( t_lower, self.t_names )

        for x in self.t_names:
            self.assertEqual( ref_alias.get_alias(x), x.lower() )
            self.assertEqual( ref_alias.get_original_name(x.lower()), x )

        with self.assertRaises(AliasError) as e:
            z = ref_alias.get_alias("non-existent")

    def test_one_to_oneness(self):

        t_len     = lambda x: str(len(x))

        with self.assertRaises(AliasError) as e:
            ref_alias = ReferenceAlias(t_len, self.t_names)

    def test_dict_input(self):
        map_dict  = {"GAPDH": "1", "VEGFA": "2", "MYC": "3"}
        ref_alias = ReferenceAlias( map_dict, self.t_names )

        self.assertEqual(ref_alias.get_alias("GAPDH"), "1")

######################################################################

def identity_alias(x):
    return x

def lower_alias(x):
    return x.lower()

def random_alias(x):
    alias_dict = {"VEGFA": "y",
                  "GAPDH": "x",
                  "MYC"  : "z"}
    return alias_dict[x]

######################################################################

class AliasTestBase(unittest.TestCase):
    def setUp(self, rename_func = None):
        self.tmp_files = list()

        self.ref_len_file       = StringIO(TRANSCRIPT_LENGTHS)
        self.annotation_file    = StringIO(TRANSCRIPT_ANNOTATION)
        self.alignment_file_1   = StringIO(READ_SET_1)
        self.alignment_file_2   = StringIO(READ_SET_2)

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

        self.sample_ribo = Ribo( self.merged_io, alias = rename_func )

    def tearDown(self):
        self.handle.close()
        self.handle_2.close()

########################################################################

class TestRiboWithNoAlias(AliasTestBase):
    """
    Should raise exceptions because no alias is defined
    """
    def test_metagene(self):
        with self.assertRaises(AliasError) as context:
            metagene_df = self.sample_ribo.get_metagene(
                             site_type       = "stop",
                             experiments     = self.sample_ribo.experiments,
                             sum_lengths     = False,
                             sum_references  = False,
                             range_lower     = 0,
                             range_upper     = 0,
                             alias           = True)

class TestRiboWithLowerAlias(AliasTestBase):

    def setUp(self):
        super().setUp(lower_alias)


    def test_metagene(self):
        metagene_df = self.sample_ribo.get_metagene(
                         site_type       = "stop",
                         experiments     = self.sample_ribo.experiments,
                         sum_lengths     = False,
                         sum_references  = False,
                         range_lower     = 0,
                         range_upper     = 0,
                         alias           = True)


        self.assertEqual( set( metagene_df.index.levels[1] ) ,
                          set( ("gapdh", "vegfa", "myc") ) )

        self.assertTrue( np.allclose (
                              metagene_df.loc["merzifon", "gapdh", 2],
                              ACTUAL_STOP_SITE_COVERAGE_length_2[0] ) )
        self.assertTrue( np.allclose (
                              metagene_df.loc["merzifon", "vegfa", 3],
                              ACTUAL_STOP_SITE_COVERAGE_length_3[1] ) )

        self.assertTrue( np.allclose (
                              metagene_df.loc["merzifon", "myc", 4],
                              ACTUAL_STOP_SITE_COVERAGE_length_4[2] ) )

    def test_region_counts(self):
        region_counts_df = self.sample_ribo.get_region_counts(
                              region_name     = "CDS",
                              sum_lengths     = True,
                              sum_references  = False,
                              range_lower     = 0,
                              range_upper     = 0,
                              experiments     = self.sample_ribo.experiments,
                              alias           = True)

        self.assertTrue( np.isclose(region_counts_df.loc["gapdh"]["merzifon"] , 13 ) )
        self.assertTrue( np.isclose(region_counts_df.loc["vegfa"]["merzifon"] , 10 ) )
        self.assertTrue( np.isclose(region_counts_df.loc["myc"]["merzifon"] , 3 ) )


    def test_coverage(self):
        coverage_df = self.sample_ribo.get_coverage(
                                        experiment  = "ankara",
                                        range_lower = 0,
                                        range_upper = 0,
                                        alias       = True)

        self.assertTrue( np.allclose(coverage_df["gapdh"],
             [0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]) )
        self.assertTrue( np.allclose(coverage_df["vegfa"],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]) )
        self.assertTrue( np.allclose(coverage_df["myc"],
             [4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]) )

    def test_coverage_transcript(self):
            gapdh_coverage_df = self.sample_ribo.get_transcript_coverage(
                                     experiment  = "ankara",
                                     transcript  = "gapdh",
                                     range_upper = 0,
                                     range_lower = 0,
                                     sum_lengths = True,
                                     alias       = True)

            self.assertTrue(np.allclose(gapdh_coverage_df.iloc[0],
                [0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ))

            myc_coverage_df = self.sample_ribo.get_transcript_coverage(
                                     experiment  = "ankara",
                                     transcript  = "myc",
                                     range_upper = 0,
                                     range_lower = 0,
                                     sum_lengths = True,
                                     alias       = True)

            self.assertTrue(np.allclose(myc_coverage_df.iloc[0],
                [4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ))




if __name__ == '__main__':
    unittest.main()
