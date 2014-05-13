# Create a consensus sequence from an alignment


import re
# read in the sequences. do this without biopython (just read in the fasta)
# collect all of column 1

infile = "/Users/sjspielman/Dropbox/Amine/HRH_Ahmad/HRH34/Round3/round3_hmm.aln"
f = open(infile, 'rU')
lines = f.readlines()
f.close()

seqList = []
for line in lines:
	find = re.search('^>', line)
	if not find:
		seqList.append(line.rstrip())

for row in range(len(seqList[0])):
	consensus = ''
	colItems = set()
	for col in range(len(seqList)):
		consensus += seqList[row][col]
		colItems.add(seqList[row][col])
	
	for item in colItems:

		
		
		
		name = lambda arguments: expression
		
		
		
		
		
		



