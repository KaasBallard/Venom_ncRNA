#!/bin/bash

# ShortStack is command line program through bioconda that according to the website: "ShortStack is a tool developed to process and analyze smallRNA-seq data with respect to a reference genome, 
# and output a comprehensive and informative annotation of all discovered small RNA genes." Basically it can discover undiscovered small RNA genes and give you annotations about these small RNAs.
# It is written in Perl, and requires samtools, RNAfold from the Vienna RNA Package in order to execute properly. When you are using it to align something to the reference genome it also requires
# gzip (if you are using zipped files, and bowtie and bowtie-build. The README for the program has a very extensive list of commands that can be used in conjuction with it. Not that it requires 4
# to 10 GB of RAM per thread, so make sure to take that into account when multithreading. It also requires quite a bit of disk space temporarily, make sure to have at least 200GB for a safe
# minimum of free space per run.
# Required inputs:
    # --genomefile GENOMEFILE : Path to the reference genome in FASTA format. Must be indexable by both samtools faidx and bowtie-build, or already indexed.
    # (--readfile [READFILE ...] | --bamfile [BAMFILE ...]) : Either --readfile or --bamfile is required.
        # --readfile [READFILE ...] : Path(s) to one or more files of reads in fastq or fasta format. May be gzip compressed. Multiple files are separated by spaces. Inputting reads triggers
        # alignments to be performed.
        # --bamfile [BAMFILE ...] : Path(s) to one or more files of aligned sRNA-seq data in BAM format. Multiple files are separated by spaces. BAM files must match the reference genome
        # given in --genomefile.

# Usage:
#ShortStack [-h] \ 
#[--version] \
#--genomefile GENOMEFILE \
#[--known_miRNAs KNOWN_MIRNAS] \
#(--readfile [READFILE ...] | --bamfile [BAMFILE ...]) \
#[--outdir OUTDIR] \
#[--adapter ADAPTER | --autotrim] \
#[--autotrim_key AUTOTRIM_KEY] \
#[--threads THREADS] \
#[--mmap {u,f,r}] \
#[--align_only] \
#[--show_secondaries] \
#[--dicermin DICERMIN [--dicermax DICERMAX] \
#[--locifile LOCIFILE | --locus LOCUS] \
#[--nohp] \
#[--dn_mirna] \
#[--strand_cutoff STRAND_CUTOFF] \
#[--mincov MINCOV] \
#[--pad PAD]

    # Required:
        #--genomfile GENOMEFILE
        #(--readfile [READFILE ...])
        # Either --readfile or --bamflie is required
    # Recommended:
        #--known_miRNAs KNOWN_MIRNAS
        #--outdir OUTDIR
        #--autotrim
        #--threads THREADS

# Code that I can reuse:

# Define variables:

    # Define genomefile path:
    Genomefile= # some path

    # Define readfile paths:
    Readfiles=(
        # I'll put something here later
    )

    # Known miRNAs:
    Known_miRNAs= # I'll put something here later

    # Output directory:
    Output_Location= # Put something here later

    # Number of threads:
    Threads=4


# This section runs the the program for each readfile:
    for ((i=0; i<${#Readfiles[@]}; i++))
    do
        # Extract output location so that the loop doesn't overwrite direcotry and creates new ones for each output instead:
        Output_Directory="${Output_Location}/$(basename "${Readfiles[$i]}")_"

        # Code for running the program:
        ShortStack --genomefile "${Genomefile}" \
            --readfile "${Readfiles[$i]}" \
            --known_miRNAs "${Known_miRNAs}" \
            --outdir "${Output_Directory}" \
            --threads "${Threads}" \
            2>&1 | tee "${Output_Directory}/output_log.txt"

        echo "Search complete for read $i"
    done