README for amine_receptors/scripts subdirectory.
Contact Stephanie at stephanie.spielman@gmail.com with questions.

## Analysis pipeline overview

1. The script run_blast.py was run using UT Austin's Center for Computational Biology and Bioinformatics cluster, Phylocluster. This script was submitted as an array job, and therefore 42 distinct file containing PSI-BLAST-collected NCBI ids are returned. The command
        ```
        cat blast_ids* >> blast_ids_all.txt
        ```
was run to place all IDs into a single file for subsequent processing.

2. These sequences were filtered using the script process_blast.py. This script performed the following tasks:
    1. Removed ID duplicates
    2. Removed non-vertebrate IDs
    3. Removed sequences with 1% ambiguities 
    4. Removed sequences annotated as either pseudogene, low-quality, and/or partial
    5. Removed sequences which GPCRHMM did not confidentally consider GPCRs.
    6. Ensures, for all protein sequences which meet criteria, correct corresponding nucleotide NCBI records exist.
    7. Obtains and saves protein and nucleotide NCBI records as well as GPCRHMM residue domain predictions.
The script produced the final output files protein_records.fasta and nucleotide_records.fasta. All corresponding in these files sequences have the same name of "protid_nucid" for convenient mapping.
This script was also run on Phylocluster, as GPCRHMM requires Linux.

4. The protein_records.fasta sequences were then aligned using the script obtain_alignment.py . 
This script performs the iterative process described in Figure 3 of the manuscript. It additionally uses the nucleotide records in nucleotide_records.fasta to convert all protein alignments to corresponding nucleotide versions.
All output files from this script, along with a description, can be found [here](../analysis/alignments/).

## Other
The script rename_tree.py was used for post-processing (i.e. figure creation) of the final phylogeny.