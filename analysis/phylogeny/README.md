README for amine_receptors/analysis/phylogeny subdirectory. 
Contact Stephanie at stephanie.spielman@gmail.com with questions.

---



## Description of Contents
_[NOTE: Phylogenetic inference and bootstrapping with performed with RAxMLv8.1.1]_

*RAxML_bestTree.inference_struc_masked_part* was built with the masked HMM-informed MSA, with two partitions for TM and EM domains.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -# 100 -s protein_struc_masked_consecpart.fasta -q partitions_consec.txt -m PROTCATLGF -n inference_struc_masked_part


*RAxML_bestTree.inference_struc_masked_nopart* was built with the masked HMM-informed MSA, with a single partition for the entire protein.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -D -p $RANDOM -# 100 -s protein_struc_masked_consecpart.fasta -m PROTCATLGF -n inference_struc_masked_nopart

       
*RAxML_bootstraps.struc_masked_part* contains bootstrap inferences for the masked HMM-informed MSA, with two partitions for TM and EM domains.
        
        raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_struc_masked_consecpart.fasta -q partitions_consec.txt -m PROTCATLGF -n struc_masked_part 

The files *RAxML_bipartitions.struc_masked_part* and *RAxML_bipartitionsBranchLabels.struc_masked_part* contain the phylogeny *RAxML_bestTree.inference_struc_masked_part* with bootstrap values incorporated.
These files were generated with the command,
```
raxmlHPC -m PROTCATLGF -f b -t RAxML_bestTree.inference_struc_masked_part -z RAxML_bootstraps.struc_masked_part -n struc_masked_part
```

__for_raxml/__    

This directory contains the specific alignment and partition files used to make trees with RAxML. These alignments' columns are *ordered differently* from the [released biogenic amine receptor alignments](../alignments/). In particular, these alignment columns have been reordered such that TM and EM partitions are consecutive. RAxML is only able to bootstrap appropriately if partitions are strictly consecutive in the alignment columns. Thus, the alignment columns needed to be reordered such that there were only 2 consecutive partitions for TM and EM, rather than alternating (following GPCR structure). *Importantly, reordering alignment columns has no bearing whatsoever on the phylogenetic analysis as all columns are treated independently.*

