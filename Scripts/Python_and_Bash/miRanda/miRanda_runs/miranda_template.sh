#/bin/bash

<<SYNOPSIS
miranda - Finds potential target sites for miRNAs in genomic sequence

Basic Usage:
miranda file1 file2 [ options ... ]
Advanced:
miranda file1 file2 [-sc score] [-en energy] [-scale scale] [-strict] [-go X] [-ge Y] [-out fileout] [-quiet] [-trim T] [-noenergy] [-restrict file]##
SYNOPSIS

<<DESCRIPTION
miRanda is an algorithm for the detection of potential microRNA target sites in genomic sequences. miRanda reads RNA sequences (such as microRNAs) from file1 and genomic DNA/RNA sequences from file2. Both of these files should be in FASTA format. This is an example of a FASTA formatted sequence:
>gi|29565487|emb|AJ550546.1| Drosophila melanogaster microRNA miR-bantam
GTGAGATCATTTTGAAAGCTG
One or more miRNA sequences from file1 are scanned against all sequences in file2 and potential target sites are reported. Potential target sites are identified using a two-step strategy. First a dynamic programming local alignment is carried out between the query miRNA sequence and the reference sequence. This alignment procedure scores based on sequence complementarity and not on sequence identity. In other words we look for A:U and G:C matches instead of A:A, G:G, etc. The G:U wobble bair is also permitted, but generally scores less than the more optimal matches. Here is an example alignment:

   Query:    3' gtCGAAAGTTTTACTAGAGTg 5' (eg. miRNA)
                  |:||||| |||||||||: 
   Ref:      5' taGTTTTCACAATGATCTCGg 3' (eg. 3'UTR)
The second phase of the algorithm takes high-scoring alignments (Those above a score threshold, defined by -sc) detected from phase 1 and estimates the thermodynamic stability of RNA duplexes based on these alignments. This second phase of the method utilizes folding routines from the RNAlib library, which is part of the ViennaRNA package written by Ivo Hofacker. At this stage we generate a constrained fictional single-stranded RNA composed of the query sequence, a linker and the reference sequence (reversed). This structure then folded using RNAlib and the minimum free energy (DG kcal/mol) is calculated for that structure.

Finally, detected targets with energies less than an energy threshold (defined by -en) are selected as potential target sites. Target site alignments passing both thresholds and other information is produced as output.
DESCRIPTION

<<OPTIONS
--help -h
Displays help, usage information and command-line options.
--version -v --license
Display version and license information.
-sc score
Set the alignment score threshold to score. Only alignments with scores >= score will be used for further analysis.
-en energy
Set the energy threshold to energy. Only alignments with energies <= energy will be used for further analysis. A negative value is required for filtering to occur.
-scale scale
Set the scaling parameter to scale. This scaling is applied to match / mismatch scores in the critical 7bp region near the 5' end of the microRNA. Many known examples of miRNA:Target duplexes are highly complementary in this region. This parameter can be thought of as a contrast function to more effectively detect alignments of this type.
-strict
Require strict alignment in the seed region (offset positions 2-8). This option prevents the detection of target sites which contain gaps or non-cannonical base pairing in this region.
-go X
Set the gap-opening penalty to X for alignments. This value must be negative.
-ge Y
Set the gap-extend penalty to Y for alignments. This value must be negative.
-out fileout
Print results to an output file called fileout.
-quiet
Quiet mode, omit notices of when scans are starting and when sequences have been loaded from input files.
-trim T
Trim reference sequences to T nucleotides. Useful when using noisy predicted 3'UTRs as reference sequences.
-noenergy
Turn off thermodynamic calculations from RNAlib. If this is used, only the alignment score threshold will be used. the -en setting will be ignored.
-restrict file
Restrict scans to those between specific miRNAs and UTRs. file should contain lines of tab separated pairs of sequence identifiers: miRNA_id <tab> target_id.
OPTIONS

<<REFERENCES
If you use this program for your research then please include the following citation:
A.J. Enright, B. John, U. Gaul, T. Tuschl, C. Sander, D.S. Marks; (2003)
MicroRNA targets in Drosophila; Genome Biology 5(1):R1.

RNAlib Citations:

I.L. Hofacker, W. Fontana, P.F. Stadler, S. Bonhoeffer, M. Tacker, P. Schuster (1994) Fast Folding and Comparison of RNA Secondary Structures. Monatshefte f. Chemie 125: 167-188

M. Zuker, P. Stiegler (1981) Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information, Nucl Acid Res 9: 133-148

J.S. McCaskill (1990) The equilibrium partition function and base pair binding probabilities for RNA secondary structures, Biopolymers 29: 1105-1119
REFERENCES

# Reuseable code:

# Define variables:

    # file1 or the RNA sequence to be compared to file2:
    RNA_sequence=(

    )

    #file2 or the genomic DNA/RNA sequence, alignment procedure uses complementarity and not sequence identity:
    genome_sequence=

    # Set alignment threshold score for usage:
    score=

    # Set energy threshold, which further restrains alignments filtered out by score
    # based on thermodynamic stability of the RNA duplexes:
    energy=

    # Set the scaling parameter:
    scale=
    
    # Gap-opening penalty, not sure what that means right now:
    X=

    # Gap-extend penalty for the alignment, not sure what that means either:
    Y=

    # This trims the reference sequences to a certain number of nucleotides:
    T=

    # This restricts the scan to certain ranges between miRNA and UTRs:
    restriction_file=

    # Name of location for the output file:
    output_location=

# Run the program:
    for ((i=0; i<${#RNA_sequence[@]}; i++))
    do
        # Extract the RNA_sequence name and create an output file name:
        output_directory="${output_location}/$(basename "${RNA_sequence[$i]}")"

        # Actually run the program:
        miranda "${RNA_sequence[$i]}" "${genome_sequence}" \
            --help \
            --version \
            -sc "${score}" \
            -en "${energy}" \
            -scale "${scale}" \
            -strict \
            -go "${X}" \
            -ge "${Y}" \
            -out "${output_directory}" \
            -quiet \
            -trim "${T}" \
            -noenergy \
            -restrict "${restriction_file}"

        echo "microRNA target site alignment complete for $i"
    done