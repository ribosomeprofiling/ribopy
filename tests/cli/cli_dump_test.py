# -*- coding: utf-8 -*-
import unittest
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
from ribopy.core.quantify import quantify_experiment

import sys
test_dir_1 = os.path.dirname(os.path.realpath(__file__))
sys.path.append(test_dir_1)
test_dir_2 = os.path.dirname(os.path.realpath(test_dir_1))
sys.path.append(test_dir_2)


from multilength_test_data import *
from ribopy.core import create_experiment

###########################################

NPROCESS = 4

ribo_output_file = "cli_create.ribo"

EXPECTED_2_and_3_COVERAGE=\
"""GAPDH\t0\t1\t1
GAPDH\t1\t2\t2
GAPDH\t3\t4\t4
GAPDH\t4\t5\t4
GAPDH\t5\t6\t1
GAPDH\t6\t7\t2
GAPDH\t7\t8\t1
GAPDH\t8\t9\t2
GAPDH\t9\t10\t1
GAPDH\t10\t11\t2
GAPDH\t13\t14\t6
GAPDH\t14\t15\t3
GAPDH\t15\t16\t3
GAPDH\t16\t17\t1
GAPDH\t18\t19\t1
GAPDH\t19\t20\t1
VEGFA\t3\t4\t5
VEGFA\t4\t5\t3
VEGFA\t5\t6\t1
VEGFA\t10\t11\t3
VEGFA\t16\t17\t3
MYC\t0\t1\t1
MYC\t6\t7\t1
MYC\t10\t11\t3"""

EXPECTED_5_COVERAGE = \
"""GAPDH\t1\t2\t1
GAPDH\t3\t4\t1
GAPDH\t5\t6\t1
GAPDH\t8\t9\t1
GAPDH\t10\t11\t1
GAPDH\t14\t15\t1
VEGFA\t3\t4\t1
VEGFA\t10\t11\t4
VEGFA\t16\t17\t2
MYC\t0\t1\t1
MYC\t10\t11\t7"""

EXPECTED_2_and_3_VEGF_TSV_COVERAGE = \
"VEGFA\t0,0,0,5,3,1,0,0,0,0,3,0,0,0,0,0,3,0,0,0,0,0"


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

        with open(REF_LEN_FILE, "w") as output_stream:
            print(TRANSCRIPT_LENGTHS, file=output_stream)
            self.tmp_files.append(REF_LEN_FILE)

        with open(ANNOTATION_FILE, "w") as output_stream:
            print(TRANSCRIPT_ANNOTATION, file=output_stream)
            self.tmp_files.append(ANNOTATION_FILE)

        with open(ALIGNMENT_FILE_1, "w") as output_stream:
            print(READ_SET_1, file= output_stream)
            self.tmp_files.append( ALIGNMENT_FILE_1 )

        create_command = ["ribopy", "create", 
                          "-a", ALIGNMENT_FILE_1,
                          "--name",          "merzifon", 
                          "--reference",     "hg38",
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
        print(output)
        self.tmp_files.append( ribo_output_file )


    def tearDown(self):
        [ os.remove(f) for f in self.tmp_files ]

#######################################################################  


class TestCLIDump(TestCLIBase):


    def test_dump_annotation(self):
        command_pieces  = ["ribopy", "dump", "annotation", ribo_output_file]
        output, error   = self.run_command(command_pieces)
        output_lines    = output.split("\n")
        expected_line_2 = ["GAPDH", "5", "15", "CDS", "0", "+"]
        output_line2    = output_lines[1].split()
        self.assertTrue( expected_line_2 == output_line2 )


    def test_dump_ref_lengths(self):
        command_pieces     = ["ribopy", "dump", "reference-lengths", 
                                ribo_output_file]
        output, error      = self.run_command(command_pieces)
        output_lines       = output.split("\n")
        expected_last_line = ["MYC", "17"]

        self.assertTrue( expected_last_line, output_lines[2] )

    def test_dump_region(self):
        #######################################
        # CDS, len 2 : 3 , no sum
        command_pieces = ["ribopy", "dump", "region-counts", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3", 
                          "--region", "CDS", 
                          ribo_output_file]

        output, error = self.run_command(command_pieces)
        output_lines  = output.split("\n")

        expected_output_lines = ["{},{},merzifon".format(DF_TRANSCRIPT, 
                                                         DF_READLENGTH),
                                 "GAPDH,2,2", "VEGFA,2,1", "MYC,2,0",
                                 "GAPDH,3,3", "VEGFA,3,2", "MYC,3,1"]

        expected_output = "\n".join(expected_output_lines)
        
        for x, y in zip( output_lines, expected_output_lines ):
            self.assertEqual(x, y)

        # CDS, len 2 : 3 , sum across transcripts
        command_pieces = ["ribopy", "dump", "region-counts", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3",
                          "--sumtrans", 
                          "--region", "CDS",  
                          ribo_output_file]

        output, error       = self.run_command(command_pieces)
        output_str          = output.strip()
        expected_output_str = "read_length,merzifon\n2,3\n3,6"

        self.assertEqual(output_str, expected_output_str)

        # CDS, len 2 : 3 , sum across lengths
        command_pieces = ["ribopy", "dump", "region-counts", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3",
                          "--sumlengths",
                          "--region", "CDS", 
                          ribo_output_file]

        output, error = self.run_command(command_pieces)
        output_str          = output.strip()
        expected_output_str = "{},merzifon\nGAPDH,5\nVEGFA,3\nMYC,1"\
                              .format(DF_TRANSCRIPT)

        self.assertEqual(output_str, expected_output_str)

        # CDS, len 2 : 3 , sum across lengths and transcripts
        command_pieces = ["ribopy", "dump", "region-counts", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3",
                          "--sumlengths",  "--sumtrans", 
                          "--region", "CDS",  
                          ribo_output_file]

        output, error = self.run_command(command_pieces)
        output_str          = output.strip()
        expected_output_str = ",merzifon\n0,9"
        self.assertEqual(output_str, expected_output_str)


    def test_dump_metagene(self):
        
        #######################################
        # CDS, len 2 : 3 , sum transcripts
        command_pieces = ["ribopy", "dump", "metagene", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3", 
                          "--site",        "start", 
                          ribo_output_file]

        output, error  = self.run_command(command_pieces)
        output_str     = output.strip() 
        expected_lines = [ ",".join( (DF_EXPERIMENT_NAME, DF_READLENGTH, 
                                        "-2", "-1", "0", "1", "2")),
                           ",".join(("merzifon", "2","1","4","1","1", "1")),
                           ",".join(("merzifon", "3","3","5","3","2","0")) ]
        expected_str = "\n".join(expected_lines)

        self.assertEqual(expected_str, output_str)    


        # CDS, len 2 : 3 , No sum
        command_pieces = ["ribopy", "dump", "metagene", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3", 
                          "--nosumtrans",
                          "--site",        "start", 
                          ribo_output_file]

        output, error  = self.run_command(command_pieces)
        output_str     = output.strip() 
        output_lines   = output_str.split()

        self.assertEqual(output_lines[0], 
                 "experiment,{},read_length,-2,-1,0,1,2"\
                       .format(DF_TRANSCRIPT))
        self.assertEqual(output_lines[6], "merzifon,MYC,3,0,0,0,0,0")


        # CDS, len 2 : 3 , Sum lengths
        command_pieces = ["ribopy", "dump", "metagene", 
                          "--experiment",  "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3", 
                          "--nosumtrans", 
                          "--sumlengths",
                          "--site", "start", 
                          ribo_output_file]

        output, error = self.run_command(command_pieces)
        output_str    = output.strip() 
        output_lines  = output_str.split()

        self.assertEqual("merzifon,VEGFA,0,5,3,1,0", output_lines[2])
        
        # CDS, len 2 : 3 , Sum lengths and transcripts on stop site

        command_pieces = ["ribopy", "dump", "metagene", 
                          "--experiment",     "merzifon",
                          "--lowerlength", "2", 
                          "--upperlength", "3", 
                          "--sumlengths",
                          "--site", "stop", 
                          ribo_output_file]

        output, error  = self.run_command(command_pieces)
        output_str     = output.strip()
        expected_lines = output_str.split()

        self.assertEqual(expected_lines[1], "merzifon,6,3,9,1,0")

        # Run without a length range
        # So this should take length min oto length max by default
        command_pieces = ["ribopy", "dump", "metagene", 
                          "--experiment", "merzifon", 
                          "--site",   "start", 
                          ribo_output_file]

        output, error  = self.run_command(command_pieces)
        output_str     = output.strip()

        expected_lines = [ ",".join( (DF_EXPERIMENT_NAME, DF_READLENGTH, 
                                        "-2", "-1", "0", "1", "2")),
                           ",".join(("merzifon", "2","1","4","1","1", "1")),
                           ",".join(("merzifon", "3","3","5","3","2","0")),
                           ",".join(("merzifon", "4","2","7","7","2","1")),
                           ",".join(("merzifon", "5","1","1","1","0","0")) ]
        expected_str   = "\n".join(expected_lines)

        self.assertEqual(expected_str, output_str) 


    def test_dump_coverage(self):
        command_pieces = ["ribopy", "dump", "coverage",
                          "--lowerlength", "2", 
                          "--upperlength", "3",
                          "--format",      "bg",
                          ribo_output_file, "merzifon" ]

        output, error  = self.run_command(command_pieces)
        output_str     = output.strip()

        self.assertEqual(output_str, EXPECTED_2_and_3_COVERAGE)

        command_pieces = ["ribopy", "dump", "coverage",
                          "--lowerlength", "5", 
                          "--upperlength", "5",
                          "--format",      "bg",
                          ribo_output_file, "merzifon" ]

        output, error  = self.run_command(command_pieces)
        output_str     = output.strip()

        self.assertEqual(output_str, EXPECTED_5_COVERAGE)





if __name__ == '__main__':
        
    unittest.main()
