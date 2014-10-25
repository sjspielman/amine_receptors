README for amine_receptors/analysis/alignments subdirectory. 
Contact Stephanie at stephanie.spielman@gmail.com with questions.

## Description of Contents
_[NOTE: all MSAs were created using protein sequences and then back-translated with the raw nucleotide sequences to create nucleotide MSAs]_

__naive__/      Contains protein and nucleotide structurally-naive MSAs, in FASTA format. 
 * __nucleotide_aln_naive.fasta__ and __protein_aln_naive.fasta__ contain nucleotide and protein MSAs, respectively.
 * __domain_alignment.fasta__ contains the MSA in protein_aln_naive.fasta but all residues have been replaced with their respective GPCR domain (O = outer, i.e. extracellular, M = membrane, I = intracellular, A = ambiguous).
 * __domain_consensus.txt__ contains simply the string showing the naive MSA's consensus domain structure.

__structural__/
 * Contains protein and nucleotide structurally-curated MSAs, in FASTA format. 
 * __nucleotide_aln_struc.fasta__ and __protein_aln_struc.fasta__ contain *unmasked* nucleotide and protein MSAs, respectively.
 * __nucleotide_aln_struc_masked.fasta__ and __protein_aln_struc_masked.fasta__ contain *masked* nucleotide and protein MSAs, respectively.
 * __domain_alignment.fasta__ contains the MSA in protein_aln_naive.fasta but all residues have been replaced with their respective GPCR domain (O = outer, i.e. extracellular, M = membrane, I = intracellular, A = ambiguous).
 * __domain_consensus.txt__ contains simply the string showing the structurally-curated MSA's consensus domain structure.