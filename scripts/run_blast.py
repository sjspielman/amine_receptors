#! /usr/bin/env python

##############################################################################
##  Script to collect biogenic amine receptor sequences using PSI-BLAST.  
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com) with some initial input from Keerthana Kumar.
##############################################################################
import sys
import re
import subprocess
from Bio import Entrez, SeqIO

# Human biogenic amine receptor NCBI IDs.
seedpairs = [ ["ADRA1A", "NP_000671.2"], ["ADRA1B", "NP_000670.1"], ["ADRA1D", "NP_000669.1"], ["ADRA2A", "NP_000672.2"], ["ADRA2B", "NP_000673.2"], ["ADRA2C", "NP_000674.2"], ["ADRB1", "NP_000675.1"], ["ADRB2", "NP_000015.1"], ["ADRB3", "NP_000016.1"], ["CHRM1", "NP_000729.2"], ["CHRM2", "NP_000730.1"], ["CHRM3", "NP_000731.1"], ["CHRM4", "NP_000732.2"], ["CHRM5", "NP_036257.1"], ["DRD1", "NP_000785.1"], ["DRD2", "NP_000786.1"], ["DRD3", "NP_000787.2"], ["DRD4", "NP_000788.2"], ["DRD1B", "NP_000789.1"], ["HRH1", "NP_000852.1"], ["HRH2", "NP_071640.1"], ["HRH3", "NP_009163.2"], ["HRH4", "NP_067637.2"], ["HTR1A", "NP_000515.2"], ["HTR1B", "NP_000854.1"], ["HTR1D", "NP_000855.1"], ["HTR1E", "NP_000856.1"], ["HTR1F", "NP_000857.1"], ["HTR2A", "NP_000612.1"], ["HTR2B", "NP_000858.2"], ["HTR2C", "NP_000859.1"], ["HTR4", "NP_001035259.1"], ["HTR5A", "NP_076917.1"], ["HTR6", "NP_000862.1"], ["HTR7", "NP_000863.1"], ["PNR", "NP_003958.2"], ["TAR1", "NP_612200.1 "], ["TAR2", "NP_001028252.1 "], ["TAR5", "NP_003958.2"], ["TAR6", "NP_778237.1"], ["TAR8", "NP_444508.1"], ["TAR9", "NP_778227.3"] ]


######### Input argument ##########



############ Blast Options ############
Entrez.email = 'stephanie.spielman@gmail.com'
e_value = 1e-20
num_iterations = 5
outfmt = '"6 sseqid pident qlen slen"'
max_target_seqs = 150000
length_error = 0.50
percent_identity = 25.0



def grab_blast_output(infile, outfile, length_error, percent_identity):
    ''' From the blast returns, create a unique set of ids, and write this set to outfile.'''
    id_set=set()
    f = open(infile, 'rU')
    for line in f:
        line = line.strip()
        acc=''
        if "ref" in line:
            pident = float(line.split()[1])
            qlen = float(line.split()[2])
            slen = float(line.split()[3])
            diff = abs(slen - qlen)/(qlen) # error computation for length differences
            if pident >= percent_identity and diff <= length_error:
                acc = line.split()[0].split("|")[3]
                if len(acc) > 5:
                    id_set.add(acc)    
    outf = open(outfile, 'w')
    for id in id_set:
        outf.write(str(id) + '\n')
    outf.close()        
    return id_set





def main():
    
    # Which seed are we running?
    usage = "Usage: python run_blast.py <rep_index>. rep_index must be in range [1,42]. "
    rep = int(sys.argv[1])
    assert(1 <= rep <= 42), usage
    
    # Prepare blastin file, which contains the seed id for the blast search.
    name = seedpairs[rep-1][0]
    seed = seedpairs[rep-1][1]
    blastin = "seed"+str(rep)+".in"
    out = open(blastin, 'w')
    out.write(seed)
    out.close()
    
    # Define output files
    blastout = "seed"+str(rep)+".out"
    idout = "blast_ids"+str(rep)+".txt"
    finalfile = "final"+str(rep)+".fasta"

    # Run the blast search.
    call_blast = "./psiblast -db nr -out "+str(blastout)+" -query "+str(blastin)+" -num_iterations "+str(num_iterations)+" -evalue "+str(e_value)+" -outfmt "+str(outfmt)
    print call_blast
    blast = subprocess.call(call_blast, shell=True)
    assert(blast == 0), "Call to BLAST failed."
    
    # Collect NCBI ids that blast returns and save to file.
    print "getting id set"
    id_set = grab_blast_output(blastout, idout, length_error, percent_identity)

main()

