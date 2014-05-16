# SJS 5/16/14.
# Simply fetch NCBI records and save to avoid future fetching. Can simultaneously cull, as desired.

import sys
import re
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "stephanie.spielman@gmail.com"

###########################################################################################################
def cullFromDesc(string):
	''' Find if has something bad.... excellent documentation, steph. really good.'''
	bad = ['LOW_QUALITY_PROTEIN', 'PSEUDOGENE', 'PARTIAL']
	string = string.upper()
	print string
	for entry in bad:
		if entry in string:
			return False
	return True
	
def getRecordAndCull(prot_id, raw_id, ncbi_dir, restrictClade = False):
	''' Retrieve NCBI record. 
		We also have the option to save only members of a clade (restrictClade).
		Also also wik, don't save anything crappy. 
	'''

	fetch = Entrez.efetch(db="protein",  id=prot_id, rettype="gb", retmode="text")
	ncbi_record = SeqIO.read(fetch, 'gb')
	keep = cullFromDesc(str(ncbi_record.description))
	if keep:
		if restrictClade:
			fullTax = ncbi_record.annotations["taxonomy"]
			if restrictClade not in fullTax:
				print "DISCARDING ", raw_id," - NOT A VERTEBRATE."
				return False
		outfile = ncbi_dir + prot_id + '.gb'
		SeqIO.write(ncbi_record, outfile, 'gb')	
	else:
		print "DISCARDING ", raw_id," - CRAPPY PROTEIN."
		return False
	return True	
###########################################################################################################


assert(len(sys.argv) == 3 or len(sys.argv) == 4), "Usage: python fetchNCBI.py <infile> <outfile> <clade restriction> . Last arg optional."

infile = sys.argv[1]
outfile = sys.argv[2]
try:	
	restrictClade = sys.argv[3]
except:
	restrictClade = False
ncbi_dir = str(os.getcwd()) + "/" + "ncbi_records/"
print ncbi_dir
if not os.path.exists(ncbi_dir):
	os.mkdir(ncbi_dir)	

records = list(SeqIO.parse(infile, 'fasta'))
outhandle = open(outfile, 'w')
for record in records:
	raw_id = str(record.id)
	seq = str(record.seq)
	find = re.search("(\wP_\d+\.\d)_(.+)$", raw_id)
	assert(find), "Bad map search."
	prot_id = find.group(1)
	
	# CHECK IF HAVE THE RECORD ALREADY!!! If the record exists locally, no need to query NCBI unnecessarily.
	if (os.path.exists(ncbi_dir + prot_id + '.gb')):
		print "have", prot_id
		outhandle.write('>'+raw_id+'\n'+seq+'\n')
	else:	
		print "dont have"
		keep = getRecordAndCull(prot_id, raw_id, ncbi_dir, restrictClade)
		if keep:
			outhandle.write('>'+raw_id+'\n'+seq+'\n')
			print "kept though!"
		else:
			print "didn't keep"
	print
outhandle.close()	
