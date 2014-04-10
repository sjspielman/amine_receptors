##### SJS, with some inspiration from KK, to run psiblast and generate fasta file with all protein sequences retrieved #####
##### Edit internal options as needed. #####

## Usage: python run_blast.py <run_index> <1 or 0> 
# run_index is a value indexed at 1 for which refseq seed id to psiblast
# the 1 or 0 is to run blast or not. If 1, run blast and then process results. if 0, just process results.
# IF RESULTS NEED TO BE PROCESSED ONLY, MUST COPY THE OUTBLAST FILE TO WDIR!!

import sys
import subprocess
from Bio import Entrez, SeqIO

def grabOutput(file, length_error, percent_identity):
	id_set=set()
	f = open(file, 'rU')
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
	return id_set

def makeFasta(id_set, finalfile, taxa):
	Entrez.email = 'stephanie.spielman@gmail.com'
	
	f = open(finalfile, "w")
	for acc in id_set:
		fetch = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=acc)
		record = SeqIO.read(fetch, "gb")
		taxonomy =  ("").join(record.annotations["taxonomy"])
		if taxa in taxonomy:
			gene = "unknown"
			feat = record.features
			for entry in feat:
				if entry.type=="Protein":
					gene = entry.qualifiers["product"][0]
					gene=re.sub(" ","_",gene)
					break
			fastatext = "> "+record.id+"_"+gene+"\n"+str(record.seq)+"\n"
			f.write(fastatext)
		else:
			continue
	f.close()		



### HUMAN AMINE RECEPTORS ###	
#seedpairs = [ ["ADRA1A", "NP_000671.2"], ["ADRA1B", "NP_000670.1"], ["ADRA1D", "NP_000669.1"], ["ADRA2A", "NP_000672.2"], ["ADRA2B", "NP_000673.2"], ["ADRA2C", "NP_000674.2"], ["ADRB1", "NP_000675.1"], ["ADRB2", "NP_000015.1"], ["ADRB3", "NP_000016.1"], ["CHRM1", "NP_000729.2"], ["CHRM2", "NP_000730.1"], ["CHRM3", "NP_000731.1"], ["CHRM4", "NP_000732.2"], ["CHRM5", "NP_036257.1"], ["DRD1", "NP_000785.1"], ["DRD2", "NP_000786.1"], ["DRD3", "NP_000787.2"], ["DRD4", "NP_000788.2"], ["DRD1B", "NP_000789.1"], ["HRH1", "NP_000852.1"], ["HRH2", "NP_071640.1"], ["HRH3", "NP_009163.2"], ["HRH4", "NP_067637.2"], ["HTR1A", "NP_000515.2"], ["HTR1B", "NP_000854.1"], ["HTR1D", "NP_000855.1"], ["HTR1E", "NP_000856.1"], ["HTR1F", "NP_000857.1"], ["HTR2A", "NP_000612.1"], ["HTR2B", "NP_000858.2"], ["HTR2C", "NP_000859.1"], ["HTR4", "NP_001035259.1"], ["HTR5A", "NP_076917.1"], ["HTR6", "NP_000862.1"], ["HTR7", "NP_000863.1"], ["PNR", "NP_003958.2"], ["TAR1", "NP_612200.1 "], ["TAR2", "NP_001028252.1 "], ["TAR5", "NP_003958.2"], ["TAR6", "NP_778237.1"], ["TAR8", "NP_444508.1"], ["TAR9", "NP_778227.3"] ]

### HUMAN HRH RECEPTORS ###
seedpairs = [ ["HRH1", "NP_000852.1"], ["HRH2", "NP_071640.1"], ["HRH3", "NP_009163.2"], ["HRH4", "NP_067637.2"] ]


############ Options ############
psiblast = "/home/sjs3495/bin/ncbi-blast-2.2.28+/bin/psiblast"
e_value = 1e-20
num_iterations = 5
outfmt = '"6 sseqid pident qlen slen"'
taxa = 'Chordata'
max_target_seqs = 10000
length_error = 0.50
percent_identity = 25.0

## Prepare blastin file. Contains the seed id.
sge = sys.argv[1] # which index we are running. be sure to minus 1 from the sge_task_id in the qsub submission script
runBlast = int(sys.argv[2])
index = int(sge) - 1 # index at 0
name = seedpairs[index][0]
seed = seedpairs[index][1]
blastin = "seed"+sge+".in"
out = open(blastin, 'w')
out.write(seed)
out.close()

## outfiles and such
blastout = "seed"+sge+".out"
finalfile = "final"+sge+".fasta"


## No matter what, results should be processed. 

if runBlast:
	# Run the blast search
	callBlast = psiblast+" -db nr -out "+str(blastout)+" -query "+str(blastin)+" -num_iterations "+str(num_iterations)+" -evalue "+str(e_value)+" -outfmt "+str(outfmt)
	print callBlast
	subprocess.call(str(callBlast), shell=True)

# Process the blast search
id_set = grabOutput(blastout, length_error, percent_identity)

# Output a fasta file with sequences
makeFasta(id_set, finalfile, taxa) # restrict to taxa














