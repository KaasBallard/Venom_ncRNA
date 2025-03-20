#/bin/bash

<<Usage
GitHub: https://github.com/ianvcaldas/smalldisco/blob/main/README.md

Smalldisco has two commands, sirna and tail, whose usage can be checked with python smalldisco.py sirna --help and python smalldisco.py tail --help, respectively.
Behind the scenes, smalldisco is actually implemented as a Snakemake pipeline and smalldisco.py is just a wrapper script. The pipeline code and helper 
scripts are in the workflow folder, and the program will not work if smalldisco.py and workflow are not in the same directory.

Overall, to run both modes of smalldisco for a full siRNA identification and tailing analysis - sirna and tail - the user needs to provide their small RNA-seq 
alignment BAM files of interest, a GTF/GFF reference annotation file, and a FASTA genome reference file. Below, we outline the details of these modes.
Usage

<<sirna_mode_overview
sirna identifies genomic regions to which user-provided reads map antisense. To this end, sirna requires:

a GTF/GFF reference annotation file
    Please be sure that your GTF/GFF file follows the standard tab-separated, nine-column format (https://useast.ensembl.org/info/website/upload/gff.html).
your small RNA-seq alignment BAM files of interest in one folder
    You can use any workflow or mapping program to create your sRNA-seq BAM files for analysis. However, be sure that your small RNA reads are trimmed of 3’ adapters (see the vignette in the "Example" section below).

The required and optional arguments can be viewed with python smalldisco.py sirna --help:
Usage: smalldisco.py sirna [OPTIONS] BAM_FOLDERs

  Find siRNA regions from antisense reads.

  This command generates putative siRNA regions based on .bam-formatted read
  alignment files in the folder BAM_FOLDERs.
sirna_mode_overview

<<sirna_Options
  -o, --out PATH                  Name of output file in BED format.
                                  [default: sirna.bed]
  -a, --annotation PATH           Genome annotation file in either GTF or GFF
                                  format.
  -k, --annotation_kind [GTF|GFF]
                                  Format of genome annotation file.  [default:
                                  GTF]
  -f, --feature TEXT              Feature type in the annotation file assumed
                                  to contain siRNA regions.  [default: CDS]
  -r X                            Minimum amount of overlapping reads to
                                  create a putative siRNA region.  [default:
                                  10; x>=1]
  -s X                            Minimum size, in base pairs, of a putative
                                  siRNA.  [default: 10; x>=1]
  --help                          Show this message and exit.

The user must define the feature type (-f) from which they wish to map antisense reads to in their GTF/GFF file. For example, if one is interested in canonical siRNAs, 
the user could define their feature type as “CDS”.

The main output file from this run will the BED formatted file sirna.bed, which will include all the identified siRNA across the samples within the BAM directory with
the source gene name as a suffix and the number of reads that support that identification of that siRNA in the last column.
sirna_Options

<<tail_mode_overview
tail uses Tailor to identify non-templated nucleotides (i.e., tails) on the 3 prime end of small RNA reads. This mode requires:

a FASTA genome reference file
your small RNA-seq alignment alignment BAM files of interest in one folder
a list of small RNA regions in BED format
The user can use the output BED file of putative siRNAs from sirna or a predefined BED file of another small RNA type.
For users interested in non-siRNA small RNAs (such as miRNAs or piRNAs): one way to create a bed file for a specific small RNA type is to obtain a GTF/GFF file for only
that small RNA type or to filter a genomic GTF/GFF file for your small RNA type of interest. Then, a tool such as gtf2bed or gff2bed 
(part of BEDOPs; https://bedops.readthedocs.io/en/latest/index.html) can be used to convert the GTF/GFF to a BED file.

NOTE: The chromosome annotation (II vs. chrII, for example) must match in the genome reference FASTA and small RNA regions BED files.
The required and optional arguments can be viewed with python smalldisco.py tail --help:

Usage: smalldisco.py tail [OPTIONS] BEDFILE BAM_FOLDERs

  Quantify tails of reads aligning to specified genome regions.

  This command quantifies tails from read alignments in .bam format found in
  BAM_FOLDERs. Only reads that overlap with certain genome regions, specified in
  .bed format in BEDFILE, are considered.
tail_mode_overview

<<tail_options
tail uses Tailor to identify non-templated nucleotides (i.e., tails) on the 3 prime end of small RNA reads. This mode requires:

a FASTA genome reference file
your small RNA-seq alignment alignment BAM files of interest in one folder
a list of small RNA regions in BED format
The user can use the output BED file of putative siRNAs from sirna or a predefined BED file of another small RNA type.
For users interested in non-siRNA small RNAs (such as miRNAs or piRNAs): one way to create a bed file for a specific small RNA type is 
to obtain a GTF/GFF file for only that small RNA type or to filter a genomic GTF/GFF file for your small RNA type of interest. Then, 
a tool such as gtf2bed or gff2bed (part of BEDOPs; https://bedops.readthedocs.io/en/latest/index.html) can be used to convert the GTF/GFF to a BED file.
NOTE: The chromosome annotation (II vs. chrII, for example) must match in the genome reference FASTA and small RNA regions BED files.

The required and optional arguments can be viewed with python smalldisco.py tail --help:
Usage: smalldisco.py tail [OPTIONS] BEDFILE BAM_FOLDERs

  Quantify tails of reads aligning to specified genome regions.

  This command quantifies tails from read alignments in .bam format found in
  BAM_FOLDERs. Only reads that overlap with certain genome regions, specified in
  .bed format in BEDFILE, are considered.

Options:
  -o, --out PATH                  Name of output file, in TSV (tab-separated
                                  values) format.  [default: tails.tsv]
  -g, --genome PATH               Reference genome in FASTA format.
  --tails-antisense / --tails-all
                                  Whether to quantify tails for antisense
                                  reads only or for all reads.  [default:
                                  tails-antisense]
  --tailor_command PATH           Path to Tailor binary executable.  [default:
                                  Tailor/bin/tailor_v1.1_linux_static]
  --tailor-min-prefix X           Minimum number of base pairs matching
                                  exactly in a read alignment before a tail
                                  can start. Equivalent to Tailor's '-l'
                                  command-line parameter.  [default: 18; x>=1]
  --help                          Show this message and exit.
The tailing result will be in the tails.tsv file, which will list the identified 3 primetail and the number of reads for each siRNA that contain that 
tail modification. If there are multiple samples within the BAM directory, the sample names will be specified in the last column.
tail_options

<<Tailor_integration
When running the tail command, we assume by default that the path to the Tailor executable is Tailor/bin/tailor_v1.1_linux_static. This will work if you are on 
Linux, are in the smalldisco repository, (e.g. cd smalldisco), and the Tailor submodule has been initialized as described in the installation section. 
If you are running smalldisco from a different folder, or have a custom Tailor installation, you must specify a path to a valid Tailor executable, for instance:
$ python smalldisco.py tail --tailor_command /usr/local/bin/tailor
Tailor_integration

# Reusable code:

# sirna mode:
# Define global variables:
  # smalldisco diretory:
  smalldisco_dir=

  # Output path:
  out_path=

  # Genome annotation file in either GTF or GFF format.
  genome_gtf=

  # Folder with all the bam files to be run
  bam_folders=

  # Minimum amount of overlapping reads to create a putative siRNA region. default 10; x>=1.
  min_overlapping_reads=

  # Minimum size in pase pairs of a putative siRNA region.
  min_size=

# Change directory to the smalldisco directory so that it runs
cd $smalldisco_dir

# Run the program:
  for files in $bam_folders
  do
    python smalldisco.py sirna \
      -o $out_path \
      -a $genome_gtf \
      -k GTF or GFF \
      -f $feature_text \
      -r $min_overlapping_reads \
      -s $min_size \
      $bam_folders \
      -help
    echo "sirna mode ran for $files"
  done

# tail mode:
# Define global variables:
  # smalldisco diretory:
  smalldisco_dir=

  # Output path:
  out_path=

  # Path to the reference genome fasta.
  genome=

  # Folder with all the bam files to be run.
  bam_folders=

  # Path to where the Tailor binary executable.
  taylor_path=

  # Minimum number of base pairs matching exactly in a read alignment before a tail can start. Equivalent to Tailor's '-l'
  tailor_min_prefix=

# Change directory to the smalldisco directory so that it runs
cd $smalldisco_dir

# Run the program:
  for file in $bam_folders
  do
    python smalldisco.py tail \
      -o PATH \
      -g $genome \
      --tails-antisense \
      --tails-all \
      --tailor_command $taylor_path \
      --tailor-min-prefix $tailor_min_prefix
      $bam_folders \
      -help
    echo "tail mode ran for $file"
  done
