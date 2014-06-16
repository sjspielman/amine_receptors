# Glue script for align, seq2gpcrhmm, consensus_structure scripts.
# GOAL: have a well-curated gpcr alignment.

import sys
import os
import subprocess
import shutil

infile = sys.argv[1] # initial inputfile, will be renamed
final_name = sys.argv[2]
hmm_dir = sys.argv[3]
dir = sys.argv[4]

sys.path.append(dir)
import gpcrhmm_fxns


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
    mafft = "linsi --quiet "+raw_file+" > "+aln_file
    runmafft = subprocess.call(mafft, shell=True)
    assert(runmafft == 0), "mafft fail"
   
    # seq2gpcrhmm
    print "seq2gpcrhmm-ing"
    gpcrhmm_fxns.amino2domain(aln_file, struc_file, hmm_dir, 0.8)

    # consensus
    print "consensus-ing"
    gpcrhmm_fxns.consStruc(struc_file, raw_file, raw_file_next, 0.1) # if 10% of a column is wrong, chuck it
    
    # check how many discarded
    print "reviewing the situation"
    num_discard = gpcrhmm_fxns.getNumDiscarded(raw_file, raw_file_next)

    print rep,"discarded",num_discard,"sequences"
    rep +=1

print "Total reps needed:", rep-1
shutil.copy(raw_file_next, final_name)
    