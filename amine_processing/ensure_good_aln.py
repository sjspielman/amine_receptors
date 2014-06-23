# Glue script for align, seq2gpcrhmm, consensus_structure scripts.
# be sure to have all pasta files in the wdir!!
# NOTE THAT SINCE STANDS ALONE OUTSIDE OF PIPELINE ON CLUSTER, IT IS DIFFERENT FROM THE HRH VERSION.
# GOAL: have a well-curated gpcr alignment.


import sys
import os
import subprocess
import shutil
import gpcrhmm_fxns

infile = sys.argv[1] # initial inputfile, will be renamed
final_name = sys.argv[2]
hmm_dir = sys.argv[3]


shutil.copy(infile, 's2g_0.fasta')
num_discard = 1 # starting value before while
rep = 0

while num_discard > 0:
    print "rep",rep
    nextrep = rep + 1
    
    raw_file = "s2g_"+str(rep)+".fasta"
    aln_file = "s2g_"+str(rep)+".aln"
    struc_file = "s2g_"+str(rep)+".struc"
    
    raw_file_next = "s2g_"+str(nextrep)+".fasta"
    
    # Align
    print "Aligning"
    shutil.copy(raw_file, "seq.fasta")
    pasta = "run_pasta.py generic_pasta_config.txt"
    runpasta = subprocess.call(pasta, shell=True)
    out_pasta = "pastajob.marker001.seq.fasta"
    shutil.copy(out_pasta, alnfile)
    
    # Remove vomit pasta files
    rm = subprocess.call("rm pastajob*", shell=True)
    
   
    # seq2gpcrhmm
    print "seq2gpcrhmm-ing"
    gpcrhmm_fxns.amino2domain(aln_file, struc_file, hmm_dir, 0.8)

    # consensus
    print "consensus-ing"
    gpcrhmm_fxns.consStruc(struc_file, raw_file, raw_file_next, 0.005) # if 0.5% of a column is wrong, chuck it
    
    # check how many discarded
    print "reviewing the situation"
    num_discard = gpcrhmm_fxns.getNumDiscarded(raw_file, raw_file_next)

    print rep,"discarded",num_discard,"sequences"
    rep +=1

print "Total reps needed:", rep-1
shutil.copy(raw_file_next, final_name)
    