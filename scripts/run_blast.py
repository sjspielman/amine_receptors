##### SJS, with some inspiration from KK, to run psiblast and save ids #####
##### Edit internal options as needed. #####

## Usage: python run_blast.py <rep>
# rep is a value indexed at 1 for which refseq seed id to psiblast

import sys
import re
import subprocess
from Bio import Entrez, SeqIO


############################## GLOBALS ###################################


######### Input argument ##########
usage = "Usage: python run_blast.py <rep_index>. For human seeds, rep_index is [1,42]. For zebrafish, it's [1,24]. "
rep = int(sys.argv[1])
assert(1 <= rep <= 42), usage


########### HUMAN AMINE RECEPTORS ##########    

# some human sequences. 42 of them.
#seedpairs = [ ["ADRA1A", "NP_000671.2"], ["ADRA1B", "NP_000670.1"], ["ADRA1D", "NP_000669.1"], ["ADRA2A", "NP_000672.2"], ["ADRA2B", "NP_000673.2"], ["ADRA2C", "NP_000674.2"], ["ADRB1", "NP_000675.1"], ["ADRB2", "NP_000015.1"], ["ADRB3", "NP_000016.1"], ["CHRM1", "NP_000729.2"], ["CHRM2", "NP_000730.1"], ["CHRM3", "NP_000731.1"], ["CHRM4", "NP_000732.2"], ["CHRM5", "NP_036257.1"], ["DRD1", "NP_000785.1"], ["DRD2", "NP_000786.1"], ["DRD3", "NP_000787.2"], ["DRD4", "NP_000788.2"], ["DRD1B", "NP_000789.1"], ["HRH1", "NP_000852.1"], ["HRH2", "NP_071640.1"], ["HRH3", "NP_009163.2"], ["HRH4", "NP_067637.2"], ["HTR1A", "NP_000515.2"], ["HTR1B", "NP_000854.1"], ["HTR1D", "NP_000855.1"], ["HTR1E", "NP_000856.1"], ["HTR1F", "NP_000857.1"], ["HTR2A", "NP_000612.1"], ["HTR2B", "NP_000858.2"], ["HTR2C", "NP_000859.1"], ["HTR4", "NP_001035259.1"], ["HTR5A", "NP_076917.1"], ["HTR6", "NP_000862.1"], ["HTR7", "NP_000863.1"], ["PNR", "NP_003958.2"], ["TAR1", "NP_612200.1 "], ["TAR2", "NP_001028252.1 "], ["TAR5", "NP_003958.2"], ["TAR6", "NP_778237.1"], ["TAR8", "NP_444508.1"], ["TAR9", "NP_778227.3"] ]

# some zebrafish sequences. 24 of them
seedpairs = [['cholin_5', 'NP_001018639.1_NM_001020803.1'], ['cholin_3', 'XP_001334664.4_XM_001334628.4'], ['cholin_4', 'XP_005169095.1_XM_005169038.1'], ['cholin_2', 'XP_005166056.1_XM_005165999.1'], ['5ht_7', 'XP_002662019.3_XM_002661973.3'], ['5ht_5A', 'NP_001007122.1_NM_001007121.2'], ['5ht_1A', 'NP_001139238.1_NM_001145766.1'], ['5ht_4', 'XP_001337671.1_XM_001337635.2'], ['5ht_6', 'XP_696681.5_XM_691589.5'], ['5ht_2B', 'NP_001038208.1_NM_001044743.1'], ['hrh_3', 'NP_001020689.1_NM_001025518.1'], ['hrh_1', 'NP_001036196.1_NM_001042731.1'], ['hrh_2', 'NP_001038803.1_NM_001045338.1'], ['adr_A-1B', 'NP_001007359.1_NM_001007358.2'], ['adr_A-2A', 'NP_997520.2_NM_207637.2'], ['trace_10', 'NP_001076507.2_NM_001083038.2'], ['trace_11', 'NP_001076546.1_NM_001083077.1'], ['trace_13E', 'NP_001076512.1_NM_001083043.1'], ['trace_12H', 'NP_001076378.1_NM_001082909.1'], ['dopa_D5', 'XP_001341592.2_XM_001341556.2'], ['dopa_D1', 'NP_001129448.1_NM_001135976.2'], ['dopa_D4', 'NP_001012638.2_NM_001012620.2'], ['dopa_D2', 'NP_922918.1_NM_197936.1'], ['dopa_D3', 'XP_005162730.1_XM_005162673.1']]

############ Blast Options ############
Entrez.email = 'stephanie.spielman@gmail.com'
e_value = 1e-20
num_iterations = 5
outfmt = '"6 sseqid pident qlen slen"'
max_target_seqs = 150000
length_error = 0.50
percent_identity = 25.0
##########################################################################



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
    ## Prepare blastin file, which contains the seed id for the blast search.
    name = seedpairs[rep-1][0]
    seed = seedpairs[rep-1][1]
    blastin = "seed"+str(rep)+".in"
    out = open(blastin, 'w')
    out.write(seed)
    out.close()

    ## Define output files
    blastout = path+"seed"+str(rep)+".out"
    idout = "blast_ids"+str(rep)+".txt"
    finalfile = path+"final"+str(rep)+".fasta"

    # Run the blast search.
    call_blast = "./psiblast -db nr -out "+str(blastout)+" -query "+str(blastin)+" -num_iterations "+str(num_iterations)+" -evalue "+str(e_value)+" -outfmt "+str(outfmt)
    blast = subprocess.call(call_blast, shell=True)
    assert(blast == 0), "Call to BLAST failed."
    
    # Collect NCBI ids that blast returns and save to file.
    print "getting id set"
    id_set = grab_blast_output(blastout, idout, length_error, percent_identity)

main()








