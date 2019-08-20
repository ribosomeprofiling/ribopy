# -*- coding: utf-8 -*-
import unittest
import os
from io import StringIO
import subprocess
import shutil
import pipes

import numpy as np
import pandas as pd
import h5py

from ribopy import create
from ribopy.core.coverage import find_coverage
from ribopy.core.get_gadgets import get_reference_names,\
                                         get_reference_lengths,\
                                         get_region_boundaries
from ribopy.settings import *
from ribopy.core.quantify import quantify_experiment
from ribopy.core import create_experiment

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)


from multilength_test_data import *
import rnaseq_data


###########################################

NPROCESS = 4

ribo_output_file = "cli_rnaseq.ribo"

RNASEQ_BED_FILE_1 = "rnaseq_1.bed"
RNASEQ_BED_FILE_2 = "rnaseq_2.bed"
RNASEQ_TSV_FILE_1 = "rnaseq_1.tsv"
RNASEQ_TSV_FILE_2 = "rnaseq_2.tsv"
RNASEQ_BAM_FILE   = "rnaseq.bam"


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
            [ (REF_LEN_FILE,      TRANSCRIPT_LENGTHS),
              (ANNOTATION_FILE,   TRANSCRIPT_ANNOTATION),
              (ALIGNMENT_FILE_1,  READ_SET_1),
              (RNASEQ_BED_FILE_1, rnaseq_data.RNASEQ_READS),
              (RNASEQ_BED_FILE_2, rnaseq_data.RNASEQ_READS_2),
              (RNASEQ_TSV_FILE_1, rnaseq_data.RNASEQ_tsv_1),
              (RNASEQ_TSV_FILE_2, rnaseq_data.RNASEQ_tsv_2) ]

        for t_file, content_str in files_and_contents:
            with open(t_file, "w") as output_stream:
                print(content_str, file = output_stream)
            self.tmp_files.append(t_file)

        if os.path.isfile(ribo_output_file):
            os.remove(ribo_output_file)

        create_command = ["ribopy", "create",
                          "--alignmentfile",  ALIGNMENT_FILE_1,  
                          "--name",           "merzifon",
                          "--reference",      "hg38",
                          "--lengths",        REF_LEN_FILE, 
                          "--annotation",     ANNOTATION_FILE,
                          "--metageneradius", METAGENE_RADIUS, 
                          "--leftspan",       LEFT_SPAN, 
                          "--rightspan",      RIGHT_SPAN,
                          "--lengthmin",      2, 
                          "--lengthmax",      5, 
                          "-n",               4,
                          ribo_output_file]

        output, error = self.run_command(create_command)
        if error:
            print("Error in creation of ribo:\n", error)
        self.tmp_files.append( ribo_output_file )


    def tearDown(self):
        [ os.remove(f) for f in self.tmp_files ]

#######################################################################  


class TestCLIRNASEQ(TestCLIBase):
    
    def test_set_rnaseq(self):
        command_pieces = ["ribopy", "rnaseq", "set",
                          "--name",   "merzifon",
                          "--alignment", RNASEQ_BED_FILE_1,
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)
        # Check that we don't get error messages        
        self.assertEqual(output, "")
        self.assertEqual(error, "")

        # Check that --force option oworks when there is
        # existing rna-seq

        command_pieces = ["ribopy", "rnaseq", "set",
                          "--name",   "merzifon",
                          "--alignment", RNASEQ_BED_FILE_2,
                          "--force",
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)
        
        self.assertEqual(output, "")
        self.assertEqual(error, "")
        
    def test_set_rnaseq_from_bam(self):
                    
        if not shutil.which('bamToBed'):
            raise FileNotFoundError("Could not find the executable bamToBed")
        if not shutil.which('bedToBam'):
            raise FileNotFoundError("Could not find the executable bedToBam")
            
        pipe_t = pipes.Template()
        pipe_t.append("bedToBam -g {}"\
                      .format( REF_LEN_FILE), "--")
      
        f = pipe_t.open(RNASEQ_BAM_FILE, 'w')
        f.write(rnaseq_data.RNASEQ_READS)
        f.close()
        
        command_pieces = ["ribopy", "rnaseq", "set",
                          "--name",      "merzifon",
                          "--alignment", RNASEQ_BAM_FILE,
                          "--format",    "bam",
                           ribo_output_file]
                           
        output, error   = self.run_command(command_pieces)

        command_pieces = ["ribopy", "rnaseq", "get",
                          "--name", "merzifon",
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)
        
        output_df = pd.read_csv(StringIO(output), 
                                index_col = (0,1),
                                header = 0,
                                sep = "\t")
        
        self.assertTrue( all(np.isclose( output_df.loc[("merzifon", "MYC"),:],
                             [0, 1, 0, 0, 2] )  ) )
        os.remove(RNASEQ_BAM_FILE)
        
    def test_set_rnaseq_from_self_table(self):
        intermediate_tsv_file = "int_rnaseq.tsv"
        
        command_pieces = ["ribopy", "rnaseq", "set",
                          "--name",   "merzifon",
                          "--alignment", RNASEQ_BED_FILE_1,
                           ribo_output_file]
        output, error  = self.run_command(command_pieces)
        
        command_pieces = ["ribopy", "rnaseq", "get",
                          "--name", "merzifon",
                          "--out", intermediate_tsv_file, 
                           ribo_output_file]                           
        output, error  = self.run_command(command_pieces)
        
        command_pieces = ["ribopy",    "rnaseq", "delete",
                          "--name",   "merzifon",
                          "--force",
                           ribo_output_file]                       
        output, error  = self.run_command(command_pieces)
        
        command_pieces = ["ribopy", "rnaseq", "set",
                          "--name",   "merzifon",
                          "--counts", intermediate_tsv_file,
                          "--force",
                           ribo_output_file]
        output, error  = self.run_command(command_pieces)
                           
        command_pieces = ["ribopy", "rnaseq", "get",
                          "--name", "merzifon",
                           ribo_output_file]
        output, error  = self.run_command(command_pieces)
        
        output_df = pd.read_csv(StringIO(output), 
                                index_col = (0,1),
                                header = 0,
                                sep = "\t")
                                
        self.assertTrue( all(np.isclose( output_df.loc[("merzifon", "GAPDH"),:],
                             [1, 0, 2, 1, 3] )  ) )
        os.remove(intermediate_tsv_file)

    #@unittest.skip("temporarily skipping tests")
    def test_get_rnaseq(self):
        command_pieces = ["ribopy", "rnaseq", "set",
                          "--name",   "merzifon",
                          "--alignment", RNASEQ_BED_FILE_1,
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)

        command_pieces = ["ribopy", "rnaseq", "get",
                          "--name", "merzifon",
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)

        output_lines          = output.split("\n")
        #VEGFA_name, VEGFA_exp = output_lines[2].split()
        ### More test here later
        output_df = pd.read_csv(StringIO(output), 
                                index_col = (0,1),
                                header = 0,
                                sep = "\t")
        
        self.assertTrue( all(np.isclose( output_df.loc[("merzifon", "VEGFA"),:],
                             [0, 1, 3, 0, 1] )  ) )
        #self.assertEqual(VEGFA_name, "VEGFA")
        #self.assertTrue(np.isclose( float(VEGFA_exp), 15.46 ))


    #@unittest.skip("temporarily skipping  tests")
    def test_delete_rnaseq(self):
        command_pieces = ["ribopy",    "rnaseq", "set",
                          "--name",   "merzifon",
                          "--alignment", RNASEQ_BED_FILE_1,
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)

        command_pieces = ["ribopy",    "rnaseq", "delete",
                          "--name",   "merzifon",
                          "--force",
                           ribo_output_file]

        output, error   = self.run_command(command_pieces)
        print(output, error)

        ribo = h5py.File(ribo_output_file)

        self.assertTrue( RNASEQ_name not in \
                         ribo[EXPERIMENTS_name]["merzifon"].keys() )

        ribo.close()


if __name__ == '__main__':
        
    unittest.main()
