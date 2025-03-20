#!/home/administrator/miniforge3/envs/biopython/bin/python

'''
Everything below here is wrong. The problem was being caused by the samtools method Yannick
showed me. Using bedtools is just better, it doesn't cause in no sequence lines and 
makes the output a single line fasta file. However, I will keep this as it was helpful for
discovering issues.

Step 2 - Convert New genomes to singline FASTAs
I noticed that the FASTA files I made are multiline, so I need to convert them to singline to avoid
any errors.

Previous step: BCFtools_get_fasta_file_three_prime_utr_2024.12.12.sh (/home/administrator/Documents/Kaas/Venom_ncRNA_project/Scripts/BCFtools/BCFtools_get_fasta_file_three_prime_utr_2024.12.12.sh)
Next step: 


'''

# Load packages
import os

# Define my conversion function
def single2multiline(fasta_directory):

    # Loop through directories and files
    for subdir, dirs, files in os.walk(fasta_directory):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.fa'):
                # Define the input and output FASTAs
                input_fasta = os.path.join(subdir, file)  # Define input FASTA path
                output_fasta = os.path.join(subdir, f"single_line_{file}")  # Define output FASTA path

                # Manually read and process the FASTA file to create single-line sequences
                with open(input_fasta, 'r') as in_file:
                    header = None
                    sequence = ''
                    headers = []
                    sequences = []
                    
                    # Read the FASTA file line by line
                    for line in in_file:
                        line = line.strip()  # Strip lines

                        # If the line is a header, i.e. it starts with '>', add it to the headers list
                        if line.startswith('>'):
                            if header:  # Save the previous header and sequence
                                headers.append(header)
                                sequences.append(sequence)
                            header = line[1:]  # Remove the '>' from the header
                            sequence = ''
                        else:
                            sequence += line  # Concatenate sequence lines

                    # Append the last header and sequence if it exists
                    if header:
                        headers.append(header)
                        sequences.append(sequence)
                
                # Write the result to the output file
                with open(output_fasta, 'w') as out_file:
                    for h, s in zip(headers, sequences):
                        out_file.write(f">{h}\n{s}\n")
                
                print(f"Single-line FASTA written to: {output_fasta}")


# Set my working directory
base_dir = '/home/administrator/Documents/Kaas/Venom_ncRNA_project/Results/bcftools/fasta_files'

# Call the function and convert FASTA files
single2multiline(fasta_directory = base_dir)