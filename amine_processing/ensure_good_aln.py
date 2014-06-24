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
if hmm_dir[-1] != '/':
    hmm_dir += '/'


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
    if not os.path.exists(aln_file):
        print "Aligning"
        mafft = "mafft "+raw_file+" > "+aln_file
        runmafft = subprocess.call(mafft, shell=True)
        assert(runmafft == 0), "mafft fail"

   
    # seq2gpcrhmm
    if not os.path.exists(struc_file):
        print "seq2gpcrhmm-ing"
        gpcrhmm_fxns.amino2domain(aln_file, struc_file, hmm_dir, 0.8)

    # consensus
    print "consensus-ing"
    gpcrhmm_fxns.consStruc(struc_file, raw_file, raw_file_next, 0.05) # if >=5% of a row is wrong, remove it

    rep +=1

print "Total reps needed:", rep-1
shutil.copy(raw_file_next, final_name)

# at 5%:
#0 discarded 145 sequences
#1 discarded 3 sequences
