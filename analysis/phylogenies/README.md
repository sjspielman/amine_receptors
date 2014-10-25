README for amine_receptors/analysis/phylogenies subdirectory. 
Contact Stephanie at stephanie.spielman@gmail.com with questions.

---



## Description of Contents
_[NOTE: All phylogenies were created with RAxMLv8.1.1]_
All phylogenies were constructed in RAxML v8.0.26 using the PROTCATLGF (LG+F matrix with CAT model of heterogeneity), and the final tree is optimized under GAMMA.
Computational resources were provided by the University of Texas at Austin's Center for Computational Biology and Bioinformatics "Phylocluster".

for_raxml/    

This directory contains the specific alignment and partition files used to make trees with RAxML. For phylogenies which used partitions (EM/TM) - RAxML is only able to bootstrap appropriately if partitions are strictly consecutive in the alignment columns. Thus, the alignment columns needed to be reordered such that there were only 2 consecutive partitions for TM and EM, rather than alternating (following GPCR structure).
*Importantly, reordering alignment columns has no bearing whatsoever on the phylogenetic analysis as all columns are treated independently.*

inferences/
    1. RAxML_bestTree.inference_aln_part       ->  raxmlHPC-PTHREADS-SSE3 -T 16 -D -p 8561 -s protein_aln_consecpart.fasta -q partitions.txt -m PROTCATLGF -n inference_aln_part 
    2. RAxML_bestTree.inference_aln_nopart     ->  raxmlHPC-PTHREADS-SSE3 -T 16 -D -p 32495 -s protein_aln_consecpart.fasta -m PROTCATLGF -n inference_aln_nopart 
    3. RAxML_bestTree.inference_masked_part    ->  raxmlHPC-PTHREADS-SSE3 -T 16 -D -p 9419 -s protein_aln_masked_consecpart.fasta -q partitions.txt -m PROTCATLGF -n inference_masked_part
    4. RAxML_bestTree.inference_masked_nopart  ->  raxmlHPC-PTHREADS-SSE3 -T 16 -D -p 15267 -s protein_aln_masked_consecpart.fasta -m PROTCATLGF -n inference_masked_nopart
    5. RAxML_bestTree.inference_naive          ->  raxmlHPC-PTHREADS-SSE3 -T 16 -D -p 27674 -s protein_aln_naive.fasta -m PROTCATLGF -n inference_naive
    
bootstraps/
    1. RAxML_bootstraps.aln_part       ->  raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_aln_consecpart.fasta -q partitions.txt -m PROTCATLGF -n aln_part
    2. RAxML_bootstraps.aln_nopart     ->  raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_aln_consecpart.fasta -m PROTCATLGF -n aln_nopart
    3. RAxML_bootstraps.masked_part    ->  raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_aln_masked_consecpart.fasta -q partitions.txt -m PROTCATLGF -n masked_part
    4. RAxML_bootstraps.masked_nopart  ->  raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_aln_masked_consecpart.fasta -m PROTCATLGF -n masked_part
    5. RAxML_bootstraps.naive          ->  raxmlHPC-PTHREADS-SSE3 -T 16 -p $RANDOM -x $RANDOM -# 200 -s protein_aln_naive.fasta -m PROTCATLGF -n naive

final_trees/
    Bootstrap values were added into phylogenies in "inferences/" by issuing this command for each:
    raxmlHPC -m PROTCATLGF -f b -t <tree_file> -z <bootstrap_file> -n <final_file>
    The resulting output file named RAxML_bipartitions.<final_file> is kept and renamed to the following- 

    1. aln_part.tre       (lnLik = -515343.618274, k_nobl = 40, k_bl = 6115, AIC = 1042917)
    2. aln_nopart.tre     (lnLik = -515991.713700, k_nobl = 20, k_bl = 6095, AIC = 1044173)
    3. masked_part.tre    (lnLik = -505500.846326, k_nobl = 40, k_bl = 6115, AIC = 1023232) ***
    4. masked_nopart.tre  (lnLik = -506397.248192, k_nobl = 20, k_bl = 6095, AIC = 1024984)
    5. naive.tre          (lnLik = -589703.683344, k_nobl = 20, k_bl = 6945, AIC = 1193279)


consec_partitions/ contains files used specifically for RAxML.

