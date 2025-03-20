<<Step8b
The point of this file is to record the commands used to generate the CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_5UTRs.gtf.
The previous step in the pipeline is: miranda_output_conversion_BED_2024-4-9.ipynb (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/miRanda/miranda_formating/miranda_output_conversion_BED_2024-4-9.ipynb)
The next steps of the pipeline is: miranda_bedtools_intersect_2024-4-9.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BEDtools/Intersect/miranda_bedtools_intersect_2024-4-9.sh)
Step8b

awk '$3 == "five_prime_utr"' /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod.gtf > /home/administrator/Documents/Kaas/Venom_ncRNA_project/Usable_data/Genome_files/CroVir_rnd1.all.maker.final.homologIDs.updatedNov2019_with_myos_geneidmod_5UTRs.gtf
