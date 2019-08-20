# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO
import subprocess

import numpy as np
import h5py

from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *
from ribopy.dump import dump_region

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)


from multilength_test_data import *

###########################################

NPROCESS = 4

ribo_output_file   = "cli_create.ribo"
ribo_output_file_2 = "merge_test.ribo"

class TestCLIBase(unittest.TestCase):

    def run_command(self, command):
      """
      Command must be a list.
      See examples below,
      """

      command_str    = tuple( map( str, command ) )
      process        = subprocess.Popen(command_str, 
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE)
      output, error  = process.communicate()
      output_str     = output.decode()
      error_str      = error.decode()

      return (output_str, error_str) 

    def setUp(self):
        self.tmp_files = list()

        files_and_contents = \
            [ (REF_LEN_FILE,     TRANSCRIPT_LENGTHS),
              (ANNOTATION_FILE,  TRANSCRIPT_ANNOTATION),
              (ALIGNMENT_FILE_1, READ_SET_1),
              (METADATA_FILE_1,  METADATA_EXPERIMENT_STR_1),
              (RIBO_META_FILE_1, RIBO_METADATA_STR_1) ]

        for t_file, content_str in files_and_contents:
            with open(t_file, "w") as output_stream:
                print(content_str, file = output_stream)
            self.tmp_files.append(t_file)

        if os.path.isfile(ribo_output_file):
            os.remove(ribo_output_file)

        create_command = ["ribopy", "create", 
                          "--alignmentfile",  ALIGNMENT_FILE_1,
                          "--name",          "merzifon", 
                          "--reference",     "hg38",
                          "--lengths",        REF_LEN_FILE, 
                          "--annotation",     ANNOTATION_FILE,
                          "--metageneradius", METAGENE_RADIUS, 
                          "--leftspan",       LEFT_SPAN, 
                          "--rightspan",      RIGHT_SPAN,
                          "--lengthmin",      2, 
                          "--lengthmax",      5, 
                          "--expmeta",        METADATA_FILE_1,
                          "--ribometa",       RIBO_META_FILE_1,
                          "-n",               4, 
                          ribo_output_file]

        output, error      = self.run_command(create_command)

        #print("Output:", output)
        self.tmp_files.append( ribo_output_file )
        if output:
            print(output)
        if error:
            print(error)


    def tearDown(self):
        [ os.remove(f) for f in self.tmp_files ]

#######################################################################        

class TestCLIMerge(TestCLIBase):
    def setUp(self):
        TestCLIBase.setUp(self)

        self.merged_ribo = "merged.ribo"

        with open(ALIGNMENT_FILE_2, "w") as output_stream:
            print(READ_SET_2, file = output_stream)
            self.tmp_files.append( ALIGNMENT_FILE_2 )


        create_command = ["ribopy", "create", 
                          "--alignmentfile",  ALIGNMENT_FILE_1,
                          "--name",           "ankara", 
                          "--reference",      "hg38",
                          "--lengths",        REF_LEN_FILE, 
                          "--annotation",     ANNOTATION_FILE,
                          "--metageneradius", METAGENE_RADIUS, 
                          "--leftspan",       LEFT_SPAN, 
                          "--rightspan",      RIGHT_SPAN,
                          "--lengthmin",      2, 
                          "--lengthmax",      5, 
                          "--expmeta",        METADATA_FILE_1,
                          "-n",               4,
                          ribo_output_file_2]

        create_command_str = list( map(str, create_command) )
        process            = subprocess.Popen(create_command_str, 
                                               stdout = subprocess.PIPE)
        output, error      = process.communicate()
        #output = output.decode()
        #print("Output:", output)

    def tearDown(self):
        TestCLIBase.tearDown(self)
        [os.remove(x) for x in (ribo_output_file_2,
                                    self.merged_ribo)]


    def test_basic_merge(self):
        merge_command = ["ribopy", "merge", self.merged_ribo, 
                         ribo_output_file, ribo_output_file_2]

        process       = subprocess.Popen(merge_command, 
                                           stdout = subprocess.PIPE)
        output, error = process.communicate()
        output        = output.decode()

        self.assertTrue( os.path.exists( self.merged_ribo ) )
        
        
    def test_attributes(self):
        merge_command = ["ribopy", "merge", self.merged_ribo, 
                         ribo_output_file, ribo_output_file_2]

        process       = subprocess.Popen(merge_command, 
                                           stdout = subprocess.PIPE)
        output, error = process.communicate()
        output        = output.decode()
        
        with h5py.File(ribo_output_file) as ribo_1, \
             h5py.File(ribo_output_file_2) as ribo_2, \
             h5py.File(self.merged_ribo) as merged_ribo:
             
             original_format_version = ribo_1.attrs[ATTRS_FORMAT_VERSION]
             merged_format_version   = merged_ribo.attrs\
                                             .get(ATTRS_FORMAT_VERSION, None)
             self.assertEqual( original_format_version, merged_format_version )
             
             original_ribo_metadata = ribo_1.attrs[USER_METADATA]
             merged_ribo_metadata   = merged_ribo.attrs\
                                             .get(USER_METADATA, None)                     
             self.assertEqual( original_ribo_metadata, merged_ribo_metadata )
             
             original_ribopy_version = ribo_1.attrs[ATTRS_VERSION]
             merged_ribopy_version   = merged_ribo.attrs\
                                             .get(ATTRS_VERSION, None)
                                             
             self.assertEqual(original_ribopy_version, merged_ribopy_version)
             
             self.assertTrue( merged_ribo.attrs.get(ATTRS_TIME, None) )
             
             CDS_counts = dump_region(ribo_handle = merged_ribo,
                             region_name = CDS_name,
                             sum_lengths = True, 
                             sum_references = False, 
                             range_lower      = 2 , 
                             range_upper      = 5 ,
                             experiment_list  = ["merzifon"])
             
             self.assertEqual(CDS_counts.loc["VEGFA", "merzifon"], 10)

        


if __name__ == '__main__':
        
    unittest.main()
