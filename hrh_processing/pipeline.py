## SJS 6/15/14. Pipeline for HRH's once sequences have been confirmed HRH's.
# Uses code in AmineReceptors/scripts/processing

import sys
import os
import re
import subprocess
from Bio import SeqIO

dir = '/Users/sjspielman/Research/AmineReceptors/hrh_processing/'
ncbi_dir = 'ncbi_records/'
hmm_dir = 'gpcrhmm_structures/'

for x in ['34']:
    infile = "confirmedGPCRs/HRH"+x+"_GPCR.fasta"

    final_prot = "final/HRH"+x+"_aa.fasta"
    final_nuc = "final/HRH"+x+"_nuc.fasta"
    final_desc = "final/HRH"+x+"_desc.fasta"


    # To begin, shorten names
    outf = open('temp0.fasta', 'w')
    records = list(SeqIO.parse(infile, 'fasta'))
    for rec in records:
        find_shortid = re.search("(\wP_\d+\.\d)_.+", str(rec.id))
        assert(find_shortid), "Couldn't parse id"
        shortid = find_shortid.group(1) 
        outf.write(">"+shortid+"\n"+str(rec.seq)+"\n")
    outf.close()
        
    # collect ncbi, if not already
    print "Collecting NCBI records"
    command1 = "python " + dir + "fetchNCBI.py temp0.fasta temp1.fasta " + ncbi_dir + " Vertebrata"
    run1 = subprocess.call(command1, shell=True)
    assert(run1 == 0)

    # remove highly similar sequences (>=98% pairwise similar)
    print "Removing highly similar sequences"
    command2 = "python " + dir + "getPairwise.py temp1.fasta temp2.fasta 0.98"
    run2 = subprocess.call(command2, shell=True)
    assert(run2 == 0)

    # obtain all gpcrhmm structures    
    print "retrieving gpcrhmm domains"
    command3 = "python " + dir + "collect_gpcrhmm.py temp2.fasta "+hmm_dir
    run3 = subprocess.call(command3, shell=True)
    assert(run3 == 0)

    # dynamic culling of the sequences whose domains won't align
    print "making sure the domains align"
    command4 = "python " + dir + "ensure_good_aln.py temp2.fasta s2g_final.fasta "+hmm_dir+" "+dir
    run4 = subprocess.call(command4, shell=True)
    assert(run4 == 0)

    # get corresponding nucleotide sequences
    print "grabbing corresponding nucleotide sequences"
    command5 = "python " + dir + "nuc_from_prot.py s2g_final.fasta " + final_prot + " " + final_nuc + " " + final_desc + " " + ncbi_dir
    run5 = subprocess.call(command5, shell=True)
    assert(run5 == 0)
    
    # Clean up    
    subprocess.call("rm temp*", shell=True)
    subprocess.call("rm aln.out", shell=True)
    subprocess.call("rm s2g*", shell=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
