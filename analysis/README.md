README for amine_receptors/analysis subdirectory.
Contact Stephanie at stephanie.spielman@gmail.com with questions.

## Description of Contents

_protein_records.fasta_ contains unaligned protein sequences input to the iterative alignment approach, outlined in Figure 1 of the paper. The file _nucleotide_records.fasta_ is its corresponding file of nucleotide sequences.
 
__gpcrhmm_records/__
 * Contains GPCR domain prediction files, generated by [GPCRHMM](http://gpcrhmm.sbc.su.se/tm.html) for all relevant proteins. Files are named using the respective NCBI protein id.

__ncbi_records/__
 * Contains NCBI genbank records for all protein and nucleotide sequences. Files are named using the respective NCBI id.

__blast/__
 * Contains two output files from the PSI-BLAST search:
  ** blast_ids_all.txt is a list of all ids returned, including duplicates
  ** blast_ids_noduplicates.txt is a list of nonduplicate ids from blast_ids_all.txt 

__alignments/__
 * Contains all sequence alignments, as well as files describing consensus structural domains in those alignments. See README within for details.

__phylogenies/__
 * Contains all phylogenetic inferences and bootstraps, as well as all RAxML-specific files. See README within for details.
