<<Step8a
The point of this file is to record the commands used to generate the CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_3UTRs_updated.gtf.
The previous step in the pipeline is: miranda_output_conversion_BED_2024-4-9.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miranda_formating/miranda_output_conversion_BED_2024-4-9.ipynb)
The next steps of the pipeline is: miranda_bedtools_intersect_2024-4-9.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BEDtools/Intersect/miranda_bedtools_intersect_2024-4-9.sh)
Step8a

# The point of this file is to record the commands used to generate the CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_3UTRs_updated.gtf so that I can fead it to BEDtools_script_for_3UTR_from_gtf_2-5-24.sh, which will give me a fasta file for all of the exons.
# Note that I have to do this becuase Aundrea/Sid made a mistake with the myotoxin gene which caused it to be skipped in the BEDtools_script_for_3UTR_from_gtf_10-25-23.sh attempt.

awk '$3 == "three_prime_utr"' /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf > /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_3UTRs_updated.gtf

<<diff_output_for_the_two_files
(bedtools) administrator@klauber:~$ diff /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_3UTRs.gtf /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_3UTRs_updated.gtf
1c1
< PE_reconstructed_10x_myo      .       three_prime_utr 1923    2076    .       +       .       gene_id "myotoxin1"; transcript_id "myotoxin_model_1"; ID "myotoxin1:three_prime_UTR"; Parent "nbis-mrna-1"; original_biotype "three_prime_UTR";
---
> PE-reconstructed-10x-myo      .       three_prime_utr 1923    2076    .       +       .       gene_id "myotoxin1"; transcript_id "myotoxin_model_1"; ID "myotoxin1:three_prime_UTR"; Parent "nbis-mrna-1"; original_biotype "three_prime_UTR";
diff
diff_output_for_the_two_files