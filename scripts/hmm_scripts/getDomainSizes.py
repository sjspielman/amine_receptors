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
		
		
		
def main(infile, outfile, hmm_dir, cutoff):
	''' Get size of each domain using gpcrhmm output.
		infile  = sequence infile
		outfile = duh?
		hmm_dir = directory to the gpcrhmm files
		cutoff  = posterior probability threshold for calling domain.
	'''
	aln = list(SeqIO.parse(infile, 'fasta'))
	outf = open(outfile, 'w')
	outf.write('name	type	full	Nterm	M1	ICL1	M2	ECL1	M3	ICL2	M4	ECL2	M5	ICL3	M6	ECL3	M7	Cterm\n')
	for record in aln:
		id = str(record.id)
		type = id.split('_')[0]
		seq = str(record.seq)
		print id
		hmmFile = hmm_dir + id + '.txt'
		total, sizeList = parseGPCRHMM(hmmFile, cutoff)
		outf.write(id + '\t' + type + '\t' + total + '\t' + '\t'.join(sizeList)+'\n')
	outf.close()
##########################################################################################	

path = os.getcwd() +'/'

if len(sys.argv) != 5:
	print "Usage: python getDomainSizes.py <infile> <outfile> <hmm_dir> <cutoff> Paths not needed, but assumes infile and hmm_dir are in same directory."
	sys.exit()

infile = path + sys.argv[1]
outfile = path + sys.argv[2]
hmm_dir = path + sys.argv[3]
if hmm_dir[-1] != '/':
	hmm_dir += '/'
cutoff = float(sys.argv[4])
main(infile, outfile, hmm_dir, cutoff)
















