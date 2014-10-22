SJS.
PIPELINE FOR AMINE RECEPTOR COLLECTION, PROCESSING, AND ALIGNMENT.

1. Run amine_receptors/scripts/run_blast.py on cluster (submit array job 1-42). This will return 42 files containing blast ids. Once all are returned, cat >> *blast_ids_all.txt* (note there were 4232 unique IDs)
2. Process the resulting blast_ids_all.txt with amine_receptors/scripts/process_blast.py. (3636 retained)
    -> Verifies NOT lowquality, pseudogene, etc.
    -> Remove sequences with >1% ambiguiuous characters.
    -> Fetch corresponding nucleotide NCBI record and confirm that it agrees with the protein sequence
    -> Save the ncbi nucleotide and protein records as well as the sequences. 
    Outfiles are found in protein_records_naive.fasta and nucleotide_records_naive.fasta. They are called "naive" because they have not been validated as GPCRs in any way whatsoever.
    All corresponding sequences have the same name of "protid_nucid" for convenient mapping.

3. Validate sequences in *_naive.fasta as GPCRs with amine_receptors/scripts/run_gpcrhmm.py . Note that this script must be run on linux with the gpcrhmm executable in the working directory!
    -> Ensures that the sequence is indeed a gpcr, using gpcrhmm. We use very conservative thresholds here.
    -> Retrieve gpcrhmm structure file from their server
    [Note: since run on cluster, produces individual fasta files for each nucleotide and protein sequence. Must cat them together.]
    Resulting files are protein_records_gpcr.fasta and nucleotide_records_gpcr.fasta. Same IDs as before. There are a total of 3464 sequences (172 "non-GPCRS" removed!)
    
4. Begin iterative alignment on protein_records_gpcr.fasta with amine_receptors/scripts/obtain_alignment.py
    -> Align.
    -> Replace sequence with domain tag. Domains called if >0.5 posterior probability, means no ambiguous (unless exactly 50-50)
    -> Remove sequences which contain >=5% sites out of line with consensus
    -> Repeat while sequences are still getting kicked out (while num_discard > 0)
    Resulting files are found in alignments/ directory (see README within for details).
