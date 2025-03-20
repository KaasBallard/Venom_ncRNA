#!/home/administrator/miniforge3/envs/biopython/bin/python

# Import packages
from Bio import SeqIO

# FASTA files to be checked
fasta_files = [
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/CDS/CV0857_viridis_North_M_genome.fasta",
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/CDS/CV0985_concolor_Other_F_genome.fasta",
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/CDS/CV0987_lutosus_Other_F_genome.fasta",
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/CDS/CV1081_viridis_Mid_M_genome.fasta",
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/CDS/CV1086_viridis_South_M_genome.fasta",
    "/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files/CDS/CV1087_viridis_North_F_genome.fasta"
]

# Define function that finds FASTA sequences of zero length
def find_zero_length_sequences(fasta_file):
    # Parse my FASTA file
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if len(record) == 0:
                print(f"Sequence {record.id} has zero length")

# Run the function in a loop
for fasta_file in fasta_files:
    find_zero_length_sequences(fasta_file)
