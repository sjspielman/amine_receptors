# Script to assess size of each domain segment using gpcrhmm results. 
# I = inner
# O = outer
# M = membrane


import re
import os
import sys
from Bio import SeqIO

def whichDomain(inner, outer, membr, cutoff):
	''' Given probabilities for each domain, which are we likely in? Lenient ambiguity, please! '''
	temp = [ inner, outer, membr ]
	which = temp.index(max(temp))
	if which == 0:
		return 'I'
	elif which == 1:
		return 'O'
	elif which == 2:
		return 'M'
	else:
		raise AssertionError("How in the hell did we get here?!!???!!!!")


##########################################################################################	
def parseGPCRHMM(file, cutoff):
	''' Parse a gpcrhmm file and return a list of domain letters to replace amino acid sequence with. ''' 
	hmm_file = open(file, 'rU')
	hmm_lines = hmm_file.readlines()
	hmm_file.close()
	
	sizeList = []
	previousDomain = 'O' # Start with O since all begin as outer/extracelluar Nterm.
	size = 0
	for line in hmm_lines:
		parseLine=re.search('^(\d+)\s\w\t(.+)\t(.+)\t(.+)\t(.+)$', line)
		if parseLine:
			total = parseLine.group(1)
			inner = float(parseLine.group(2))
			outer =  float(parseLine.group(3)) + float(parseLine.group(5))
			membr = float(parseLine.group(4))
			currentDomain = whichDomain(inner, outer, membr, cutoff)
			if currentDomain == previousDomain:
				size += 1
			else:
				sizeList.append(str(size))
				previousDomain = currentDomain
				size = 1
	sizeList.append(str(size))
	return (total, sizeList)
		
		
		
def main(infile, hmm_dir, cutoff):
	''' Get size of each domain using gpcrhmm output.
		infile  = sequence infile
		hmm_dir = directory to the gpcrhmm files
		cutoff  = posterior probability threshold for calling domain.
	'''
	aln = list(SeqIO.parse(infile, 'fasta'))
	print 'name	full	Nterm	M1	ICL1	M2	ECL1	M3	ICL2	M4	ECL2	M5	ICL3	M6	ECL3	M7	Cterm'
	for record in aln:
		id = str(record.id)
		seq = str(record.seq)
		hmmFile = hmm_dir + id + '.txt'
		total, sizeList = parseGPCRHMM(hmmFile, cutoff)
		print id + '\t' + total + '\t' + '\t'.join(sizeList)
##########################################################################################	


infile = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH2/HRH2_map.fasta"
hmm_dir = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH2/gpcrhmm_HRH2_map/"
cutoff = 0.5
main(infile, hmm_dir, cutoff)
















