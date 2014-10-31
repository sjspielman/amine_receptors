README for amine_receptors/analysis/phylogenies subdirectory. 
Contact Stephanie at stephanie.spielman@gmail.com with questions.

---



## Description of Contents
_[NOTE: All phylogenies were created with RAxMLv8.1.1]_


__inferences/__ 

Contains phylogeny inferences (without bootstrap values!), built with the commands given below.

1. *RAxML_bestTree.inference_naive* was built with the structurally-naive MSA.
       
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -s protein_struc_naive.fasta -m PROTCATLGF -n inference_naive
2. *RAxML_bestTree.inference_struc_nopart* was built with the unmasked structurally-curated MSA, with a single partition.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -s protein_struc_consecpart.fasta -m PROTCATLGF -n inference_struc_nopart 
3. *RAxML_bestTree.inference_struc_part* was built with the unmasked structurally-curated MSA, with two partitions for TM and EM domains.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -s protein_struc_consecpart.fasta -q partitions.txt -m PROTCATLGF -n inference_struc_part 
4. *RAxML_bestTree.inference_masked_nopart* was built with the masked structurally-curated MSA, with a single partition.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -s protein_struc_masked_consecpart.fasta -m PROTCATLGF -n inference_masked_nopart
5. *RAxML_bestTree.inference_masked_part* was built with the masked structurally-curated MSA, with two partitions for TM and EM domains.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -s protein_struc_masked_consecpart.fasta -q partitions.txt -m PROTCATLGF -n inference_masked_part


__inferences/__ 

Contains bootstrap trees, built with the commands given below.

1. *RAxML_bootstraps.naive* was built with the structurally-naive MSA.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_struc_naive.fasta -m PROTCATLGF -n naive
2. *RAxML_bootstraps.struc_nopart* was built with the unmasked structurally-curated MSA, with a single partition.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_struc_consecpart.fasta -m PROTCATLGF -n struc_nopart 
3. *RAxML_bootstraps.struc_part* was built with the unmasked structurally-curated MSA, with two partitions for TM and EM domains.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_struc_consecpart.fasta -m PROTCATLGF -n struc_part 
4. *RAxML_bootstraps.masked_nopart* was built with the masked structurally-curated MSA, with a single partition.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_struc_masked_consecpart.fasta -m PROTCATLGF -n masked_nopart 
5. *RAxML_bootstraps.masked_part* was built with the masked structurally-curated MSA, with two partitions for TM and EM domains.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_struc_masked_consecpart.fasta -m PROTCATLGF -n masked_part 


__final_trees/__

Contains final phylogenetic inferences. Bootstraps were added onto each inference using the general command, 
```
raxmlHPC -m PROTCATLGF -f b -t <tree_file> -z <bootstrap_file> -n <final_file> 
```
The resulting output file named RAxML_bipartitions.<final_file> was retained and renamed to the following - 

1. struc_part.tre
2. struc_nopart.tre
3. masked_part.tre
4. masked_nopart.tre
5. naive.tre


__for_raxml/__    

This directory contains the specific alignment and partition files used to make trees with RAxML. These alignments' columns are *ordered differently* from the [released biogenic amine receptor alignments](../alignments/). In particular, these alignment columns have been reordered such that TM and EM partitions are consecutive. RAxML is only able to bootstrap appropriately if partitions are strictly consecutive in the alignment columns. Thus, the alignment columns needed to be reordered such that there were only 2 consecutive partitions for TM and EM, rather than alternating (following GPCR structure). *Importantly, reordering alignment columns has no bearing whatsoever on the phylogenetic analysis as all columns are treated independently.*

