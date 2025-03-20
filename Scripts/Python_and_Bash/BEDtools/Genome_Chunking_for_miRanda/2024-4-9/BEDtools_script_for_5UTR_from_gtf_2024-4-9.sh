
<<Step4b
The point of this file is to record the commands used to generate the Cvv_2017_genome_with_myos_5UTR.fasta that I then fed to miranda.
The previous steps of this pipeline is: mature_miRNA_extraction_2024-4-9.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/ShortStack_out_to_miRanda_data/mature_miRNA_extraction_2024-4-9.sh)
The next steps of this pipeline is: miranda_mature_2024-4-9.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miRanda_runs/2024-4-9/miranda_mature_2024-4-9.sh).
Step4b

# I used this to get the Cvv_2017_genome_with_myos_5UTR.fasta file.
bedtools getfasta -fo /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos_5UTR.fasta -fi /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/Cvv_2017_genome_with_myos.fasta -bed /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_5UTRs.gtf
