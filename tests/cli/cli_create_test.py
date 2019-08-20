# -*- coding: utf-8 -*-
import unittest
from unittest import mock
from unittest.mock import patch
import os
from io import StringIO
import subprocess
import shutil

import numpy as np
import h5py
import pipes

from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)


from multilength_test_data import *

###########################################

NPROCESS = 4

ribo_output_file = "cli_create.ribo"


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

        # We needed this adjustment for bedToBam compatiblitiy
        # Our bed file has whitespace instead of tabs
        # this makes bedToBam complain
        self.TRANSCRIPT_LENGTHS_str = TRANSCRIPT_LENGTHS.replace(" ","\t")
        self.READ_SET_1_str         = READ_SET_1.replace("   ", "\t")
        self.READ_SET_1_str         = self.READ_SET_1_str.replace("  ", "\t")
        self.READ_SET_1_str         = self.READ_SET_1_str.replace(" ", "\t")


        files_and_contents = \
            [ (REF_LEN_FILE,     self.TRANSCRIPT_LENGTHS_str),
              (ANNOTATION_FILE,  TRANSCRIPT_ANNOTATION),
              (ALIGNMENT_FILE_1, self.READ_SET_1_str),
              (METADATA_FILE_1,  METADATA_EXPERIMENT_STR_1),
              (RIBO_META_FILE_1, RIBO_METADATA_STR_1) ]

        for t_file, content_str in files_and_contents:
            with open(t_file, "w") as output_stream:
                print(content_str, file = output_stream)
            self.tmp_files.append(t_file)

        if os.path.isfile(ribo_output_file):
            os.remove(ribo_output_file)

        create_command = ["ribopy", "create", 
                          "--name",           "merzifon",
                          "--alignmentfile",  ALIGNMENT_FILE_1, 
                          "--reference",      "hg38",
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
        pass
        [ os.remove(f) for f in self.tmp_files ]

#######################################################################        

#@unittest.skip("temporarily skipping plot tests")
class TestCLICreate(TestCLIBase):

    ######################################################
    def test_create_ribo_exists(self):
        self.assertTrue(os.path.exists(ribo_output_file))
        myribo             = h5py.File(ribo_output_file, "r")
        experiments_handle = myribo[EXPERIMENTS_name]
        experiment_names   = tuple(experiments_handle.keys())
        myribo.close()

        self.assertTrue("merzifon" in experiment_names)



    def test_info_display(self):
        command_pieces = [ "ribopy", "info" , ribo_output_file ]
        output, error  = self.run_command(command_pieces)
        if error:
            print("error: ", error)

        self.assertTrue( "merzifon" in output )
        self.assertTrue( "Version" in output )


    def test_exp_name_from_path(self):
        """
        If the name is not provided,
        it can be extracted from the filepath. 
        """

        ribo_experiment_file = "this_file.ribo"
        self.tmp_files.append(ribo_experiment_file)

        create_command = ["ribopy", "create", 
                          "--alignmentfile",  ALIGNMENT_FILE_1, 
                          "--reference",      "hg38",
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
                          ribo_experiment_file]

        output, error  = self.run_command(create_command)

        with h5py.File(ribo_experiment_file , "r") as ribo:
            experiment_list = tuple( ribo[EXPERIMENTS_name].keys() )
            self.assertTrue( "this_file" in experiment_list )


    def test_create_from_bam_file(self):
        """
        Test if we can create ribo files directly from bam files
        """
        if not shutil.which('bamToBed'):
            raise FileNotFoundError("Could not find the executable bamToBed")
        if not shutil.which('bedToBam'):
            raise FileNotFoundError("Could not find the executable bedToBam")

        ribo_experiment_file = "from_bam_file.ribo"
        self.tmp_files.append(ribo_experiment_file)

        ALIGNMENT_FILE_BAM_1 = "frombed.bam"
        self.tmp_files.append(ALIGNMENT_FILE_BAM_1)

        # First we make the bam file
        pipe_t = pipes.Template()
        pipe_t.append("bedToBam -g {}"\
                      .format( REF_LEN_FILE), "--")
      
        f = pipe_t.open(ALIGNMENT_FILE_BAM_1, 'w')
        f.write(self.READ_SET_1_str)
        f.close()
        
        if not os.path.isfile(ALIGNMENT_FILE_BAM_1):
            raise FileNotFoundError("Couldnt convert bed file to bam.")

        create_command = ["ribopy", "create", 
                          "--name",            "merzifon",
                          "--alignmentfile",   ALIGNMENT_FILE_BAM_1,
                          "--alignmentformat", "bam", 
                          "--reference",       "hg38",
                          "--lengths",         REF_LEN_FILE, 
                          "--annotation",      ANNOTATION_FILE,
                          "--metageneradius",  METAGENE_RADIUS, 
                          "--leftspan",        LEFT_SPAN, 
                          "--rightspan",       RIGHT_SPAN,
                          "--lengthmin",       2, 
                          "--lengthmax",       5, 
                          "--expmeta",         METADATA_FILE_1,
                          "--ribometa",        RIBO_META_FILE_1,
                          "-n",                4,
                          ribo_experiment_file]

        output, error  = self.run_command(create_command)

        with h5py.File(ribo_experiment_file , "r") as ribo:
            experiment_list = tuple( ribo[EXPERIMENTS_name].keys() )
            self.assertTrue( "merzifon" in experiment_list )



class TestCLIMetadata(TestCLIBase):

    
    #######################################################


    def test_ribo_metadata(self):
        command_pieces = [ "ribopy", "metadata" , "get", 
                            ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        self.assertTrue("bowtie2" in output)


    def test_delete_ribo_metadata(self):
        command_pieces = [ "ribopy", "metadata" , "delete", 
                            "--force", ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        command_pieces = [ "ribopy", "metadata" , "get", 
                            ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        self.assertTrue("bowtie2" not in output)


    def test_experiment_metadata(self):
        command_pieces = [ "ribopy", "metadata" , "get", 
                           "--name", "merzifon",
                           ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        self.assertTrue("HeLa" in output)


    def test_delete_experiment_metadata(self):
        command_pieces = [ "ribopy", "metadata" , "delete", 
                           "--name", "merzifon",
                            "--force", ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        command_pieces = [ "ribopy", "metadata" , "get", 
                           "--name", "merzifon",
                            ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        self.assertTrue("HeLa" not in output)
        self.assertTrue(not output)
        self.assertTrue(not error)


    def test_change_ribo_metadata(self):

        new_ribo_meta_file = "new_ribo_meta_file.yaml"
        with open(new_ribo_meta_file, "w") as meta_out_stream:
            print("new_aligner: bwa", file = meta_out_stream)

        self.tmp_files.append(new_ribo_meta_file)

        command_pieces = [ "ribopy", "metadata" , "set",
                           "--meta",   new_ribo_meta_file, 
                            "--force", ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        command_pieces = [ "ribopy", "metadata" , "get", 
                            ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        self.assertTrue("bwa" in output)

    def test_change_experiment_metadata(self):

        new_exp_meta_file = "new_exp_meta_file.yaml"
        with open(new_exp_meta_file, "w") as meta_out_stream:
            print("new_type: HEK293", file = meta_out_stream)

        self.tmp_files.append(new_exp_meta_file)

        command_pieces = [ "ribopy", "metadata" , "set",
                           "--meta", new_exp_meta_file, 
                           "--name", "merzifon",
                            "--force", ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        command_pieces = [ "ribopy",  "metadata" , "get", 
                           "--name", "merzifon",
                            ribo_output_file ]
        output, error  = self.run_command(command_pieces)

        self.assertTrue("HEK293" in output)



###########################################################

if __name__ == '__main__':
        
    unittest.main()
