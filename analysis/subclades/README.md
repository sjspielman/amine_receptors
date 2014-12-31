README for amine_receptors/subclades subdirectory.
Contact Stephanie at stephanie.spielman@gmail.com with questions.

---


## Description of Contents

1. <subclade>_nuc.fasta and <subclade>_aa.fasta are nucleotide and protein, respectively, alignments from the maked structurally-informed alignment.
These alignments contain only those receptors in the respective <subclade>. Columns which are only gaps in these subclade alignments have been removed.

2. <subclade>.tree contains the corresponding subclade from the structurally-partitioned ML phylogeny.

3. <subclade>.domains contains a single string giving each residue's domain (outer, inner, membrane). Note that, as gap-only columns were removed from each sub-alignment, the domain positions differ accordingly from the full alignment's domain consensus. 