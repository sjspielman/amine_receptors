amine_receptors
==============
This repository accompanies the paper [*Comprehensive, structurally-informed alignment and phylogeny of vertebrate biogenic amine receptors*](https://peerj.com/preprints/571/), by Spielman SJ, Kumar K, and Wilke CO. 
Contact Stephanie at stephanie.spielman@gmail.com with any questions.

## Quick links!
1. Structurally-informed, *unmasked* [protein MSA](./analysis/alignments/structural/protein_aln_struc.fasta) and [nucleotide MSA](./analysis/alignments/structural/nucleotide_aln_struc.fasta), both in FASTA format.
 
2. Structurally-informed, *masked* [protein MSA](./analysis/alignments/structural/protein_aln_struc_masked.fasta) and [nucleotide MSA](./analysis/alignments/structural/nucleotide_aln_struc_masked.fasta), both in FASTA format.

3. [ML partitioned-phylogeny](./analysis/phylogeny/RAxML_bipartitions.struc_masked_part) created with structurally-informed *masked* alignment, in newick format with branch lengths and bootstrap support included.

4. [Description](./analysis/sequence_descriptions.txt) of all sequences in final structurally-informed alignment. This is a tab-delimited file containing the ProteinID/NucleotideID, gene, and taxonomy.

5. [Alignments and phylogenies](./analysis/subclades/) for each primary biogenic amine receptors subclade. Created from the masked structurally-informed MSA and its corresponding partitioned ML phylogeny.

## Description of Contents 
__analysis/__
 * Contains data (incl. sequences, alignments, phylogenies, sequence descriptions) used and generated during analysis pipeline. See README within for details.

__scripts/__
 * Contains all code used during analysis pipeline. See README within for details.

__Manuscript/__
 * Contains manuscript, figures, tables.

