#!/bin/bash

# sirna mode:
# Define global variables:
  # smalldisco diretory:
  smalldisco_dir="/home/administrator/smalldisco"

  # Output path:
  out_path="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/smalldisco/2024-4-27/sirna"

  # Make director if it doesn't exist
  mkdir -p "$out_path"

  # Genome annotation file in either GTF or GFF format.
  genome_gtf="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf"

  # Folder with all the bam files to be run
  bam_folders=(
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Trimmed_sRNA_BAMs"
  )

  # Minimum amount of overlapping reads to create a putative siRNA region. default 10; x>=1.
  #min_overlapping_reads=

  # Minimum size in pase pairs of a putative siRNA region.
  #min_size=

# Change directory to the smalldisco directory so that it runs
cd $smalldisco_dir

# Run the program:
  for files in $bam_folders
  do
    python smalldisco.py sirna -o $out_path -a $genome_gtf $bam_folders 
    echo "sirna mode ran for $files"
  done