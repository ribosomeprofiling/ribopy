# -*- coding: utf-8 -*-

from .main import *
from ..create import create_ribo_file_wrapper

###########################################################

@cli.command()
@click.argument('ribofile',      type = click.Path( ))
@click.option('-a', '--alignmentfile',
              help     = "Aligned Ribo-Seq Data File\n"
                         "If no file is provided, "
                         "it is read from standard input."
                         "Alignment file can be in bed or bam format."
                         "Bed files with '.gz' or '.gzip'" 
                         "extension are assumed to be gzipped." 
                         "So they will be automatiocally unzipped."
                         "If no aligment file is given," 
                         "it is read from the standard input."
                         "Thus, the user can pipe the output of some process"
                         "to ribopy to produce ribo files. ",
              type     = click.Path(exists = True),
              required = False)
@click.option('-f', '--alignmentformat' ,
               help         = "Alignment File Format",
               type         = click.Choice(['bed', 'bam']),
               default      = "bed",
               show_default = True,
               required     = False)
@click.option('--name',
               help     = "Experiment name\n"
                          "If name is not provided, experiment "
                          "name is extracted from the ribo file path. "
                          "If --name parameter is omitted, experiment name is determined "
                          "using the ribo file path. For example if the ribo file path is "
                          "/home/user/data/WT.ribo, "
                          "then, experiment name is set to 'WT'. Experiment name can only contain "
                          "alphanumeric characters, '_', '-' and '.'. ",
               type     = click.STRING ,
               required = False)
@click.option('--reference',
               help     = "Reference name\n"
                          "It is good practice to give a unique "
                          "name to a particular "
                          "transcriptome and annotation pair and "
                          "use it consistently."
                          "For convenience, the user has the freedom" 
                          "to choose a name for the transcript reference "
                          "and the annotation, used to create a ribo file,"
                          "If used consistently, one can tell whether "
                          "two ribo files are coming from the same annotation or not."
                          "Note that ribo files having different reference names "
                          "can not be merged.",
               type     = click.STRING ,
               required = True)
@click.option('--lengthsfile', '--lengths',
              help     = "A tab-separated file containing"
                         " ref. name and ref. lengths",
              type     = click.Path(exists = True),
              required = True)
@click.option('--annotationfile', '--annotation',
              help     = "A bed file defining UTR5,"
                         " CDS and UTR3 regions."
                         "Each transcript must be annotated. "
                         "More explicitly, the coordinates of the regions UTR5, " 
                         "CDS and UTR3 must be provided for each transcript. "
                         "This annotation is provided in a BED file. "
                         "The annotation can not contain gaps. "
                         "All nucleotides of a transcript must belong to a region. "
                         "If a transcript does not have a CDS or UTR3 region, "
                         "it must be excluded from the lengths file and annotation. "
                         "Note that bed files are 0-based, "
                         "start position is included "
                         "and the end position is excluded. "
                         "The order and transcript names in the lengths file "
                         "must match the order in the annotation file. "
                         "Annotation data is kept in ribo files. "
                         "It is possible to extract the annotation from a ribo file "
                         "using the dump command.",
              type     = click.Path(exists = True),
              required = True)
@click.option('--metageneradius',  '--radius',
              help     = "Number of nucleotides on either side of start"
                         " / stop sites for metagene analysis.\n\n"
                         "Aligning transcripts by their start / stop sites"
                         "and aggregating the coverage gives us metagene data."
                         "The number of nucleotides to be taken on either side of"
                         "start / stop site is given by the metagene radius.",
              type     = click.INT,
              default  = 50,
              required = False)
@click.option('--leftspan', '-l',
              help         = "Number of nucleotides to the left of start"
                             " / stop sites for defining UTR5 / UTR3"
                             " junction regions.",
              type         = click.INT,
              default      = 35,
              show_default = True,
              required     = False)
@click.option('--rightspan', '-r',
              help         = "Number of nucleotides to the right of start"
                             " / stop sites for defining "
                             "UTR5 / UTR3 junction regions.",
              type         = click.INT,
              default      = 15,
              show_default = True,
              required     = False)
@click.option('--lengthmin', '--min',
              help         = "Minimum read length to be counted",
              type         = click.INT,
              default      = 15,
              show_default = True,
              required     = False)
@click.option('--lengthmax', '--max',
              help         = "Maximum read length to be counted",
              type         = click.INT,
              default      = 35,
              show_default = True,
              required     = False)
@click.option('--ribometa',
              help     = "Metadata file, for the ribo file, in yaml format.",
              type     = click.Path(exists = True),
              required = False)
@click.option('--expmeta',
              help     = "Metadata file, for the experiment, in yaml format.",
              type     = click.Path(exists = True),
              required = False)
@click.option( '--nocoverage',
               is_flag = True,
               help    = "Do not store coverage. "
                         "By default, coverage IS stored "
                         "at nucleotide resolution, "
                         "of each transcript, for each read length. "
                         "Turning this flag on decreases size of ribo file "
                         "at the cost of coverage data. "
                         "From the alignment data, five prime end of the reads, "
                         "mapping to each nucleotide position, "
                         "of each transcript is computed for each read length. "
                         "We call this coverage. "
                         "Coverage is used for metagene analysis and region counts. "
                         "By default, ribopy stores coverage data in the ribo file. "
                         "Keeping coverage increases the ribo file size considerably. "
                         "Users can choose NOT to keep coverage by setting this flag. ")
@click.option('--nprocess', '-n',
              help         = "Number of cores to be used.",
              type         = click.INT,
              default      = 1,
              show_default = True,
              required     = False)

def create(ribofile,
           alignmentfile,
           alignmentformat,
           lengthsfile,
           name,
           reference,
           annotationfile,
           metageneradius,
           leftspan,      rightspan,
           lengthmin,     lengthmax,
           nocoverage,
           ribometa,
           expmeta,
           nprocess):
    """Creates a ribo file from a given reference, annotation and alignment file.

    The resulting ribo file contains a single experiment
    whose name is provided in "--name". If "--name" is not provided,
    the experiment name will be extracted from the ribo file path
    after stripping the file extension.

    RiboPy works on transcriptomic coordinates.
    In other words, all coordinates are relative to transcripts
    where the first nucleotide of the transcript is always zero.
    Each entry in the alignment reference must come from a single transcript.
    Ribo-Seq data aligned against genomic coordinates is NOT usable
    with RiboPy.

    A note on transcript <-> reference:
     Sequencing reads are mapped against a reference
     to generate alignment files in bed or bam format. 
     Each entry of this reference is coming
     from a transcript. Therefore, in this context, refrence and transcript
     correspond to the same entity. Thus, the names 'reference' and 'transcript'
     are used interchangebly. For example 'reference names' and
     'transcript names' refer to the same list of names.

    Reference Lengths File:
     Reference lengths must be provided in the 'lengths' option.
     This must be a tab separated file where transcript names
     are in the first coulmn and
     the transcript lengths are in the seocond column.

         \b
      
         Example:
          ============    ====
          TRANSCRIPT_1    1512
          TRANSCRIPT_2    1387
          ============    ====


    Annotation Examples:
     Below are some valid and invalid annotation file examples.
     ribopy will report an error and fail to generate a ribo file
     if an invalid annotation file is given.

         \b
     
         VALID ANNOTATION:
          (Assuming length of TRANSCRIPT_1 is 1252)
           
           ============  ====  ====  =====   ==  ===
           TRANSCRIPT_1  0     21     UTR5    0   \+
           TRANSCRIPT_1  21    1041   CDS     0   \+
           TRANSCRIPT_1  1041  1252   UTR3    0   \+
           ============  ====  ====  =====   ==  ===

         \b
          
         INVALID ANNOTATION:
          (Assuming length of TRANSCRIPT_2 is 1000.
          Nucleotide positions 32,33 and 34 are not annotated)
           
           ============  ===  =====   =====  ==  ===
           TRANSCRIPT_1  0     32     UTR5    0  \+
           TRANSCRIPT_1  35    920    CDS     0  \+
           TRANSCRIPT_1  920   1000   UTR3    0  \+
           ============  ===  =====   =====  ==  ===

         \b
         
         INVALID ANNOTATION:
          (Assuming length of TRANSCRIPT_3 is 1200)
          (There is no CDS region.)
     
           ============  ==  ======   =====  ==  ==
           TRANSCRIPT_3  0     50     UTR5    0  \+
           TRANSCRIPT_3  50    1200   UTR3    0  \+
           ============  ==  ======   =====  ==  ==

         \b
         
         INVALID ANNOTATION:
          (Assuming length of TRANSCRIPT_4 is 1000)
          (There is no UTR3 region.)
     
           ============  ===  ======   =====  ==  ==
           TRANSCRIPT_3  0     700     UTR5    0  \+
           TRANSCRIPT_3  100    1000   CDS     0  \+
           ============  ===  ======   =====  ==  ==

    Region Counts & Left / Right Span:
     For each transcript, the number of reads aligning to
     each region (UTR5, CDS, and UTR3) are stored.
     For this quantification, we exclude nucleotides in some
     proximity of start and stop sites.

     This proximity is defined by taking
     'leftspan' many nucleotides to the left of start / stop sites and
     'rightspan' many nucleotides to the right of start / stop sites.
     The regions around start and stop sites are called
     UTR5_junction and UTR3_junction respectively.
      \b
      
       UTR5_juction: Nucleotides 'around' start site.
                     This region is between UTR5 and CDS.
      \b
      
       UTR3_juction: Nucleotides 'around' stop site.
                     This region is between CDS and UTR3.
                     
       Note that 'around' is precisely defined by 'leftspan' and 'rightspan'
       arguments.


    Read Length Range:
     Quantified Ribo-Seq data (metagene coverage, region counts,
     transcript ocoverage (if any))
     is stored for each read length for a given range.
     This range is defined by --lengthmin and --lengthmax.
     Both values are inclusive.
     So, for example, for --lengthmin 19 --lengthmax 21,
     region counts are computed and stored for
     RNA fragments of length 19,20 and 21 separetely.
     For human Ribo-Seq data, a range from 15 to 35
     can be sufficient for conventional experiments.

    Metadata:
     Metadata, either for the ribo file or for any experiment,
     is an optional argument. Metadata must be provided in yaml format.

     \b
     Users can provide metadata for
        i)  ribo file
           --ribometa
        ii) experiment
           --expmeta
    
     Metadata is provided in pairs of the form 'label: value'

        \b
                
        Example Metadata:
         ============  ======
          cell-line:   HEK
          enzyme:      RNASEI
         ============  ======


    Examples:
     Below are some ribo file generation examples in different settings.
     
     1) Some command line exaples to create ribo files.
     Create a file named WT.ribo in the current directory.
     Provide alignment data in a zipped bed file.
     Do NOT store coverage data.

            \b 
            
            .. code:: bash      
             
             ribopy create --name WT \\
                --alignmentfile WT.bed.gz \\
                --reference appris_human_v1 \\
                --lengths appris_len.tsv \\
                --annotation appris_regions.bed \\
                --radius 50 \\
                -l 35 -r 15 \\
                --lengthmin 15 --lengthmax 35 \\
                --nocoverage \\
                WT.ribo     
             
             
 

     2) Create a file named Treatment_1.ribo in the current directory.
     Provide alignment data in a zipped bam file.
     Store coverage data.

     \b
     
     .. code:: bash
     
             ribopy create --name Treatment_1 \\
                 --alignmentfile Treatment_1.bam \\
                 --format bam \\
                 --reference appris_human_v1 \\
                 --lengths appris_len.tsv \\
                 --annotation appris_regions.bed \\
                 --radius 50 \\
                 -l 35 -r 15 \\
                 --lengthmin 15 --lengthmax 35 \\
                Treatment_1.ribo
                
     3) Create a file named WT_2.ribo in the current directory.
     Provide alignment data in a zipped bed file.
     Do NOT store coverage data.
     Read metadata of this experiment from WT_2_meta.yaml

     \b
                 
     .. code:: bash 
     
             ribopy create --name WT_2 \\
                 --alignmentfile WT_2.bed.gz \\
                 --reference appris_human_v1 \\
                 --lengths appris_len.tsv \\
                 --annotation appris_regions.bed \\
                 --radius 50 \\
                 -l 35 -r 15 \\
                 --lengthmin 15 --lengthmax 35 \\
                 --nocoverage \\
                 --expmeta WT_2_meta.yaml \\
                 WT_2.ribo  
    """

    print("creating the ribo file {}...".format(ribofile))

    create_ribo_file_wrapper(
                     ribo_file_path      = ribofile,
                     experiment_name     = name,
                     alignment_file      = alignmentfile,
                     alignment_format    = alignmentformat,
                     reference_name      = reference,
                     lengths_file        = lengthsfile,
                     annotation_file     = annotationfile,
                     metagene_radius     = metageneradius,
                     left_span           = leftspan,
                     right_span          = rightspan,
                     length_min          = lengthmin,
                     length_max          = lengthmax,
                     ribo_metadata       = ribometa,
                     experiment_metadata = expmeta,
                     store_coverage      = not nocoverage,
                     nprocess            = nprocess)

    print("Done.")
