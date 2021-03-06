README for amine_receptors/analysis subdirectory.
Contact Stephanie at stephanie.spielman@gmail.com with questions.

## Description of Contents

__protein_records.fasta__ contains unaligned protein sequences input to the iterative alignment approach, outlined in [Figure 3](../Manuscript/figures/alignment_flowchart.pdf) of the paper. The file __nucleotide_records.fasta__ is its corresponding file of nucleotide sequences. All corresponding sequences have the same id: \<protein.id_nucleotide.id\>, for example XP_004175317.1_XM_004175269.1 .
 
__gpcrhmm_records/__
 * Contains GPCR domain prediction files, generated by [GPCRHMM](http://gpcrhmm.sbc.su.se/tm.html), for all relevant protein sequences. Files are named using the respective NCBI protein id.

__ncbi_records/__
 * Contains NCBI genbank records for all protein and nucleotide sequences. Files are named using the respective NCBI id.

__blast/__
 * Contains two output files from the PSI-BLAST search:
     * blast_ids_all.txt is a list of all ids returned, including duplicates
     * blast_ids_noduplicates.txt is a list of nonduplicate ids from blast_ids_all.txt 

__alignments/__
 * Contains all sequence alignments, as well as files describing consensus structural domains in those alignments. See README within for details.

__phylogeny/__
 * Contains phylogenetic inference (and bootstraps) for the structurally-informed masked protein alignment. See README within for details.

