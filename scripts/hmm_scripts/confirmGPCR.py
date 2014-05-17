## SJS 5/16/14.
## Script to check if sequences are all really GPCRs, using gpcrhmm.
## MUST BE RUN ON LINUX TO USE GPCRHMM!!

import sys
import subprocess
import os
import re
from Bio import SeqIO


assert( len(sys.argv) == 4), "Usage: python confirmGPCR.py <infile> <outfile> <outfile_directory>. infile should be in working directory."

gpcrhmm_path = 'gpcrhmm/'
wdir = os.getcwd() + '/'
infile = wdir + sys.argv[1] # file containing all sequences. should be fasta format
outfile = sys.argv[2]
rdir = sys.argv[3]
outfile = rdir + '/' + outfile


###################################################################################
def runGPCRHMM(gpcrhmm_path, infile, outfile):
	command = 'perl '+gpcrhmm_path+'gpcrhmm.pl '+infile+' > '+outfile
	run = subprocess.call(command, shell=True)
	assert(run == 0), "gpcrhmm did not properly run."

def parseGPCRHMM(seqid, outfile):
	''' Second line contains the information '''
	outf = open(outfile, 'rU')
	line = outf.readlines()[1]
	outf.close()
	getScores = re.search('^'+seqid+'\s+(-*\d*\.*\d*)\s+(-*\d*\.*\d*)\s+(\w+)\s*', line)
	assert(getScores), "Couldn't parse the gpcrhmm output file."
	globalScore = getScores.group(1)
	localScore = getScores.group(2)
	pred = getScores.group(3)
	if pred == 'GPCR':
		if globalScore == '-':
			globalScore = 0.
		if localScore == '-':
			localScore = 0.
		return(True, float(globalScore), float(localScore))
	else:
		return (False, globalScore, localScore)
###################################################################################

globalThresh = 10
localThresh = 10

gpcrhmmIn = wdir+'tempin.fasta'
gpcrhmmOut = wdir+'tempout.txt'

records = list(SeqIO.parse(infile, 'fasta'))
outf = open(outfile, 'w')
outculled = open('removed.txt', 'w')
outculled.write('id\tglobalScore\tlocalScore\n')
for record in records:
	seq = str(record.seq)
	id = str(record.id)
	print id
	temp = open(gpcrhmmIn, 'w')
	temp.write(">"+id+"\n"+seq+"\n")
	temp.close()
	runGPCRHMM(gpcrhmm_path, gpcrhmmIn, gpcrhmmOut)
	isGPCR, globalScore, localScore = parseGPCRHMM(id, gpcrhmmOut)
	if isGPCR:
		print "it's a gpcr with",globalScore, localScore
		if globalScore >= globalThresh and localScore >= localThresh:
			print "keep!"
			outf.write(">"+id+"\n"+seq+"\n")
		else:
			print "but discard"
			outculled.write(id+'\t'+str(globalScore)+'\t'+str(localScore)+'\n')
	else:
		print "not a gpcr at all"
		outculled.write(id+'\t'+str(globalScore)+'\t'+str(localScore)+'\n')
	print
outf.close()
















	
	
	
	
	
	
	
	
	
	
