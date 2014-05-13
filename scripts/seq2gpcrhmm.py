# Script to replace sequence with domain
# I = inner
# O = outer
# M = membrane
# A = ambiguous (<0.7)

import re
import os
import sys
from Bio import AlignIO

##########################################################################################	
def probs2letter(inner, outer, membr, cutoff):
	''' Turn probabilities going down the protein into a letter signifying domain.
		I = inner. O = outer. M = membrane. A = ambiguous.
		Return list of letters to replace alignment amino acids with.
	'''
	domain = []
	for i in range(len(inner)):
		temp = [ inner[i], outer[i], membr[i] ]
		if max(temp) >= cutoff:		
			which = temp.index(max(temp))
			if which == 0:
				domain.append('I')
			elif which == 1:
				domain.append('O')
			elif which == 2:
				domain.append('M')
			else:
				raise AssertionError("How in the hell did we get here?!!???!!!!")
		else:
			domain.append('A')
	assert(len(domain) == len(inner)), "A domain was not identified for all alignment positions. Womp womp."
	return domain			

def parseGPCRHMM(file, cutoff):
	''' Parse a gpcrhmm file and return a list of domain letters to replace amino acid sequence with. ''' 
	hmm_file = open(file, 'rU')
	hmm_lines = hmm_file.readlines()
	hmm_file.close()
	
	innerList = []
	outerList = []
	membrList = []
	for line in hmm_lines:
		parseLine=re.search('^\d+\s\w\t(.+)\t(.+)\t(.+)\t(.+)$', line)
		if parseLine:
			innerList.append( float(parseLine.group(1)) )
			outerList.append( float(parseLine.group(2)) + float(parseLine.group(4)) )
			membrList.append( float(parseLine.group(3)) )
	assert(len(innerList) == len(outerList) == len(membrList) ), "GPCRHMM file not properly parsed."
	newLetters = probs2letter(innerList, outerList, membrList, cutoff)
	return newLetters	

def replaceSeq(aaseq, letters):
	''' Replace the amino acid string (which includes gaps!!) with the domain letters. '''
	nogap_count = 0
	newseq = ''
	for position in aaseq:
		if position == '-':
			newseq += '-'
		else:
			newseq += letters[nogap_count]
			nogap_count += 1
	assert(len(newseq) == len(aaseq)), "AA sequence not properly replaced with domain indicators."
	return newseq

def amino2domain(infile, outfile, hmm_dir, cutoff):
	''' Replace a GPCR alignment amino acids with domain.
		infile  = alignment infile
		outfile = outfile with new domain letters in lieu of amino acids
		hmm_dir = directory to the gpcrhmm files
		cutoff  = posterior probability threshold for calling domain.
	'''
	aln = AlignIO.read(infile, 'fasta')
	outaln = open(outfile, 'w')
	for record in aln:
		id = str(record.id)
		seq = str(record.seq)
		hmmFile = hmm_dir + id + '.txt'
		print id
		newLetters = parseGPCRHMM(hmmFile, cutoff)
		finalSeq = replaceSeq(seq, newLetters)
		outaln.write(">"+id+"\n"+finalSeq+"\n")
	outaln.close()
##########################################################################################	


infile = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/Round3/SubAligns/HRH4_mammals.aln"
outfile = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/Round3/SubAligns/HRH4_mammals_hmm.aln"
hmm_dir = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/Round3/gpcrhmm_round3/"
cutoff = 0.8
amino2domain(infile, outfile, hmm_dir, cutoff)

### WONKY H3 MAMMALS ###
## >H3_Talpidae.Condylura1 is totally fucked up.
## >H3_Tenrecinae.Echinops1
##>H3_Caprinae.Capra1
## >H3_Bovinae.Bos2
## >H3_Marmotini.Ictidomys1
## >H3_Camelidae.Vicugna1, >H3_Camelidae.Camelus1
## >H3_Cricetinae.Cricetulus1
## >H3_Equidae.Equus1


## The HRH4 mammals are simply gorgeous. Except for the gibbon?? >H4_Hylobatidae.Nomascus1 which is apparently a low-quality protein, so we should probably discard it anyways!!!














































gpcrhmm_dir = ''