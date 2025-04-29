#!/bin/bash

<<Step0.4-Finding_BPP
The point of this file is to record how I found where BPP is in the ISOSeq data so that I can get it's CDS, 3'UTR, and 5'UTR.
The previous step of this pipeline is: seqtk_2024-5-13.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/seqtk/seqtk_2024-5-13.sh)
The next step of this pipeline is: Venom_BPP_Best_Blast_Hits_2024-5-13.R (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/blast/Venom_BPP_Best_Blast_Hits_2024-5-13.R)
Step0.4-Finding_BPP

# Define variables
    # IsoSeq data
    isoseq="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/IsoSeq/Cvv_VG_CLR_hg_transcripts.fasta"

    # Venom BPP sequence
    bpp="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Venom_BPP/BPP_transcript_model.fasta"

    # Create path for the database
    isoseq_database="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/IsoSeq/Cvv_VG_CLR_hg_transcripts"

    # Create out file path
    outfile="/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/blast/Venom_BPP/Venom_BPP_Blast_Hits_2024-5-13.txt"

    # Cut off for e-values
    evalue="1e-60"

# Run makeblastdb to create IsoSeq database
makeblastdb -in "$isoseq" -dbtype nucl -out "$isoseq_database"

# Run blastn with the following options
blastn -task megablast -query "${bpp}" -db "${isoseq_database}" -evalue "${evalue}" -out "${outfile}" -outfmt 6
