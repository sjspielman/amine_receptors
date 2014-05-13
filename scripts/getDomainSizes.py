# Script to assess size of each domain segment using gpcrhmm results. 
# I = inner
# O = outer
# M = membrane


import re
import os
import sys
from Bio import AlignIO

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
	domainList = []
	previousDomain = 'O' # Start with O since all begin as outer/extracelluar Nterm.
	size = 0
	for line in hmm_lines:
		parseLine=re.search('^(\d+)\s\w\t(.+)\t(.+)\t(.+)\t(.+)$', line)
		if parseLine:
			#print "LINE:", parseLine.group(1)
			inner = float(parseLine.group(2))
			outer =  float(parseLine.group(3)) + float(parseLine.group(5))
			membr = float(parseLine.group(4))
			currentDomain = whichDomain(inner, outer, membr, cutoff)
			if currentDomain == previousDomain:
				size += 1
			else:
				#print "switch", previousDomain 
				domainList.append(previousDomain)
				sizeList.append(str(size))
				previousDomain = currentDomain
				size = 1
	domainList.append(previousDomain)
	sizeList.append(str(size))
	print '\t'.join(domainList)
	print '\t'.join(sizeList)
	#assert 1==0
	#return sizeList
		
		
		
def main(infile, outfile, hmm_dir, cutoff):
	''' Replace a GPCR alignment amino acids with domain.
		infile  = alignment infile
		outfile = outfile with new domain letters in lieu of amino acids
		hmm_dir = directory to the gpcrhmm files
		cutoff  = posterior probability threshold for calling domain.
	'''
	aln = AlignIO.read(infile, 'fasta')
	#outaln = open(outfile, 'w')
	for record in aln:
		id = str(record.id)
		seq = str(record.seq)
		hmmFile = hmm_dir + id + '.txt'
		print id
		parseGPCRHMM(hmmFile, cutoff)
		print
		#finalSeq = replaceSeq(seq, newLetters)
		#outaln.write(">"+id+"\n"+finalSeq+"\n")
	#outaln.close()
##########################################################################################	


infile = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/ROUND1/HRH34.aln"
outfile = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/ROUND1/HRH34_hmm.aln"
hmm_dir = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/ROUND1/gpcrhmm_HRH34/"
cutoff = 0.5
main(infile, outfile, hmm_dir, cutoff)
















