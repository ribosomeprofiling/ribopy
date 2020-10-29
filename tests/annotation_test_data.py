# ARTIFICIAL test data for testing ribopy


REF_LEN_FILE = "ref_lengths.tsv"
TEST_RIBO = "test_create.ribo"
ANNOTATION_FILE = "annotation.bed"
ALIGNMENT_FILE_1 = "alignment_file_1.bed"

###########################################
METAGENE_RADIUS = 2
LEFT_SPAN       = 3
RIGHT_SPAN      = 2
NPROCESS        = 1
LENGTH_MIN      = 2
LENGTH_MAX      = 6

SINGLE_TRANSCRIPT_LENGTHS=\
"""Gene_1   20"""


SINGLE_TRANSCRIPT_ANNOTATION=\
"""Gene_1    0    6    UTR5    0    +
Gene_1    6    15    CDS    0    +
Gene_1    15    20    UTR3    0    +"""

# This should fail
SINGLE_TRANSCRIPT_ANNOTATION_2=\
"""Gene_1    0    6    UTR5    0    +
Gene_1    6    5    CDS    0    +
Gene_1    5    20    UTR3    0    +"""

NOUTR5_SINGLE_TRANSCRIPT_ANNOTATION=\
"""Gene_1    0    9    CDS    0    +
Gene_1    9    20    UTR3    0    +"""

NOUTR3_SINGLE_TRANSCRIPT_ANNOTATION=\
"""Gene_1    0    6    UTR5    0    +
Gene_1    6    20    CDS    0    +"""

#################################################

GENERIC_TRANSCRIPT_LENGTHS = \
"""Gene_1   20
Gene_NO-UTR5    18
Gene_NO-UTR3    17
Gene_Only-CDS   15
Gene_short_UTR5 12
Gene_short_UTR3 22"""

GENERIC_TRANSCRIPT_ANNOTATION = \
"""Gene_1    0    6    UTR5    0    +
Gene_1    6    15    CDS    0    +
Gene_1    15    20    UTR3    0    +
Gene_NO-UTR5    0    12    CDS    0    +
Gene_NO-UTR5    12    18    UTR3    0    +
Gene_NO-UTR3    0    7    UTR5    0    +
Gene_NO-UTR3    7    17    CDS    0    +
Gene_Only-CDS    0    15    CDS    0    +
Gene_short_UTR5    0    1    UTR5    0    +
Gene_short_UTR5    1    8    CDS    0    +
Gene_short_UTR5    8    12    UTR3    0    +
Gene_short_UTR3    0    8    UTR5    0    +
Gene_short_UTR3    8    20    CDS    0    +
Gene_short_UTR3    21    22    UTR3    0    +"""

#################################################

READ_SET_1=\
"""Gene_NO-UTR5 0  5  Gene_NO-UTR5_read_1 0   +
Gene_NO-UTR5 1  5  Gene_NO-UTR5_read_2 0    +
Gene_NO-UTR5 4  8  Gene_NO-UTR5_read_3 0    +
Gene_NO-UTR5 15  18  Gene_NO-UTR5_read_4 0    +
Gene_NO-UTR3 0  4  Gene_NO-UTR3_read_1 0   +
Gene_NO-UTR3 3  7  Gene_NO-UTR3_read_2 0   +
Gene_NO-UTR3 6  10  Gene_NO-UTR3_read_3 0   +
Gene_NO-UTR3 7  11  Gene_NO-UTR3_read_4 0   +
Gene_NO-UTR3 8  10  Gene_NO-UTR3_read_5 0   +
Gene_NO-UTR3 11  14  Gene_NO-UTR3_read_6 0   +
Gene_NO-UTR3 12  15  Gene_NO-UTR3_read_7 0   +
Gene_NO-UTR3 15  17  Gene_NO-UTR3_read_8 0   +
Gene_Only-CDS 0  5  read_1 0   +
Gene_Only-CDS 0  4  read_2 0   +
Gene_Only-CDS 1  6  read_3 0   +
Gene_Only-CDS 2  6  read_4 0   +
Gene_Only-CDS 3  6  read_5 0   +
Gene_Only-CDS 4  6  read_6 0   +
Gene_Only-CDS 5  6  read_7 0   +
Gene_Only-CDS 10  16  read_8 0   +
Gene_Only-CDS 10  12  read_9 0   +
Gene_Only-CDS 11  17  read_10 0   +
Gene_Only-CDS 13  17  read_11 0   +
Gene_short_UTR3 19  21  read_1 0   +"""

READ_SET_2=\
"""GAPDH 250  276  read_1 0   +
BRCA 1  33  read_2 0 +"""


GENERAL_METADATA=\
{
    "metagene_radius": 5,
    "left_span": 3,
    "right_span": 2
}
