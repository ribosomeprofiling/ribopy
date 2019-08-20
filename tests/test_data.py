# ARTIFICIAL test data for testing ribopy


REF_LEN_FILE = "ref_lengths.tsv"
TEST_RIBO = "test_create.ribo"
ANNOTATION_FILE = "annotation.bed"
ALIGNMENT_FILE_1 = "alignment_file_1.bed"

###########################################


TRANSCRIPT_LENGTHS=\
"""GAPDH	1290
VEGFA	875
MYC	1462
BRCA	565"""


TRANSCRIPT_ANNOTATION=\
"""GAPDH    0   50  UTR5    0   -
GAPDH   50  1130    CDS 0   -
GAPDH   1130    1290 UTR3   0   -
VEGFA   0   675 CDS 0   +
VEGFA   675 875 UTR3    0   +
MYC 0   90  UTR5    0   -
MYC 90  1251    CDS 0   -
MYC 1251    1462    UTR3    0   -
BRCA    0   2   UTR5    0   +
BRCA    2   561 CDS 0   +
BRCA    561 565 UTR3    0   +"""

READ_SET_1=\
"""GAPDH 250  276  read_1 0   +
BRCA 1  33  read_2 0 +"""

GENERAL_METADATA=\
{
    "metagene_radius": 5,
    "left_span": 3,
    "right_span": 2
}

RNASEQ_DATA_1 = \
"GAPDH\t2.89\nVEGFA\t15.46\nMYC\t8"
