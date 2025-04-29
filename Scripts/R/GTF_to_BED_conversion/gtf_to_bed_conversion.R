# Last Edited: 2024/12/11

# Step 1: Converting GTF to BED
# The point of this script is to convert the genome .gtf into a bed file for all of the three_prime_utr, five_prime_utr, and CDS sequences.
# Next Step: BCFtools_get_fasta_file_three_prime_utr_2024.12.12.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BCFtools/BCFtools_get_fasta_file_three_prime_utr_2024.12.12.sh)

#### Set up ####

# Load packages
library(tidyverse)

# Set working directory
setwd('/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files')

# Set variable for gtf
gtf_file <- 'CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP.gtf'

# Load the gtf into memory
gtf <- read.table(
    file = gtf_file, 
    sep = '\t', 
    header = FALSE,
    quote = '' # Prevent stripping of quotes so the regex can work
) %>% 
    rename(
        seqid = V1, source = V2, type = V3, start = V4, end = V5, score = V6,
        strand = V7, phase = V8, attributes = V9
    )
glimpse(gtf)


#### Make 3'UTR .bed file ####

# Get a three prime utr data frame
gtf_3utr <- gtf %>% 
    filter(type == 'three_prime_utr') %>% 
    # Change the base offset to format as a bed file
    mutate(
        start = start - 1,
        transcript_id = str_extract(attributes, 'transcript_id "[^"]+"') %>% 
                        str_remove_all('transcript_id "') %>% 
                        str_remove_all('"')  # Extract and clean transcript_id
    ) %>% 
    # Select the columns of a bed file
    # Select and rename columns for BED file format
    select(
        chrom = seqid, chromStart = start, chromEnd = end, name = transcript_id, score, strand
    )
View(gtf_3utr)

# Write to file
write.table(
    gtf_3utr,
    file = 'CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_three_prime_utr.bed',
    quote = FALSE,
    sep = '\t',
    row.names = FALSE,
    col.names = FALSE
)

#### Make 5'UTR .bed file ####

# Get a three prime utr data frame
gtf_5utr <- gtf %>% 
    filter(type == 'five_prime_utr') %>% 
    # Change the base offset to format as a bed file
    mutate(
        start = start - 1,
        transcript_id = str_extract(attributes, 'transcript_id "[^"]+"') %>% 
                        str_remove_all('transcript_id "') %>% 
                        str_remove_all('"')  # Extract and clean transcript_id
    ) %>% 
    # Select the columns of a bed file
    # Select and rename columns for BED file format
    select(
        chrom = seqid, chromStart = start, chromEnd = end, name = transcript_id, score, strand
    )
View(gtf_5utr)

# Write to file
write.table(
    gtf_5utr,
    file = 'CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_five_prime_utr.bed',
    quote = FALSE,
    sep = '\t',
    row.names = FALSE,
    col.names = FALSE
)


#### Make CDS .bed file ####

# Get a three prime utr data frame
gtf_cds <- gtf %>% 
    filter(type == 'CDS') %>% 
    # Change the base offset to format as a bed file
    mutate(
        start = start - 1,
        transcript_id = str_extract(attributes, 'transcript_id "[^"]+"') %>% 
                        str_remove_all('transcript_id "') %>% 
                        str_remove_all('"')  # Extract and clean transcript_id
    ) %>% 
    # Select the columns of a bed file
    # Select and rename columns for BED file format
    select(
        chrom = seqid, chromStart = start, chromEnd = end, name = transcript_id, score, strand
    )
View(gtf_cds)

# Write to file
write.table(
    gtf_cds,
    file = 'CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_edited_with_BPP_CDS.bed',
    quote = FALSE,
    sep = '\t',
    row.names = FALSE,
    col.names = FALSE
)
