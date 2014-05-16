## SJS 5/16/14.
## Script to check if sequences are all really GPCRs, using gpcrhmm.
## MUST BE RUN ON LINUX TO USE GPCRHMM!!

import sys
import subprocess
import os
import re
from Bio import SeqIO


assert( len(sys.argv) == 3), "Usage: python confirmGPCR.py <infile> <outfile>. File should be in working directory. No paths needed."

gpcrhmm_path = '/home/sjs3495/bin/gpcrhmm/'
wdir = os.getcwd() + '/'
infile = wdir + sys.argv[1] # file containing all sequences. should be fasta format
outfile = wdir + sys.argv[2]

###################################################################################
def runGPCRHMM(gpcrhmm_path, infile, outfile):
	command = 'perl '+gpcrhmm_path+'/gpcrhmm.pl '+infile+' > '+outfile
	run = subprocess.call(command, shell=True)
	assert(run == 0), "gpcrhmm did not properly run."

def parseGPCRHMM(seqid, outfile):
	''' Second line contains the information '''
	outf = open(outfile, 'rU')
	line = outf.readlines()[1]
	outf.close()
	getScores = re.search('^'+seqid+'\s+(\d*\.*\d*)\s+(\d*\.*\d*)\s+\w+', line)
	assert(getScores), "Couldn't parse the gpcrhmm output file."
	globalScore = float(getScores.group(1))
	localScore = float(getScores.group(2))
	return (globalScore, localScore)
###################################################################################

globalThresh = 50
localThresh = 50

gpcrhmmIn = wdir+'tempin.fasta'
gpcrhmmOut = wdir+'tempout.txt'

records = list(SeqIO.parse(infile, 'fasta'))
outf = open(outfile, 'w')

for record in records:
	seq = str(record.seq)
	id = str(record.id)
	print id
	temp = open(gpcrhmmIn, 'w')
	temp.write(">"+id+"\n"+seq+"\n")
	temp.close()
	runGPCRHMM(gpcrhmm_path, gpcrhmmIn, gpcrhmmOut)
	globalScore, localScore = parseGPCRHMM(id, gpcrhmmOut)
	
	if globalScore >= globalThresh and localScore >= localThresh:
		outf.write(">"+id+"\n"+seq+"\n")
	else:
		print id, globalScore, localScore
outf.close()
















	
	
	
	
	
	
	
	
	
	