# SJS 5/16/14.
# Map so names are now taxonomy and receptor type.

import sys
import re
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "stephanie.spielman@gmail.com"

# Broad classification.
amine_receptors = {"DOPAMINE": "DOPA", "HISTAMINE": "HRH", "TRACE": "TRACE", "SEROTONIN": "5HT", "5-HYDROXYTRYPTAMINE": "5HT", "ADRENERGIC": "ADREN", "ADRENOCEPTOR": "ADREN", "MUSCARINIC": "CHOLIN", "ACETYLCHOLINE": "CHOLIN"}

###########################################################################################################
def getTaxonomy(ncbi_dir, prot_id):
	ncbi_record = SeqIO.read(ncbi_dir+prot_id+'.gb', 'gb')
	feat = ncbi_record.features
	taxonomy = ''
	for entry in feat:
		if entry.type == 'source':
			temp = entry.qualifiers['organism'][0]
			new = temp.split(' ')
			for entry in new:
				camel = entry[0].upper() + entry[1:].lower()
				taxonomy += camel
				taxonomy += '_'
	taxonomy = taxonomy.rstrip('_')
	assert(taxonomy != ''), "Couldn't get taxonomy'"
	return taxonomy

###########################################################################################################



###########################################################################################################
def newMap(infile, outfile, mapfile, ncbi_dir, dict):
	
	records = list(SeqIO.parse(infile, 'fasta'))
	outseqs = open(outfile, 'w')
	map = open(mapfile, 'w')
	
	have_dict = {}
	for record in records:
		original = str(record.id)
		find = re.search("(\wP_\d+\.\d)_(.+)$", original)
		assert(find), "Bad map search."
		protid = find.group(1)
		desc = find.group(2).upper()
		found=False
		rtype = 'UNKNOWN'
		for entry in dict:
			if entry in desc:
				rtype = dict[entry]
		taxa = getTaxonomy(ncbi_dir, protid)
			
		newid = rtype+'_'+taxa
		if newid not in have_dict:
			have_dict[newid] = 1
			newid += '1'
		else:
			have_dict[newid] += 1
			newid = newid + str(have_dict[newid])
		print newid
		outseqs.write(">"+newid+"\n"+str(record.seq)+"\n")
		map.write(newid+'\t'+protid+'\n')

	outseqs.close()
	map.close()
###########################################################################################################


if len(sys.argv) != 4:
	print
	print "Usage: python <infile> <outfile> <mapfile>" 
	print
	print "NOTE: infile contains input sequences, outfile will contain sequences w/ new names, and mapfile will be a tdf to map back as needed."
	sys.exit()

infile = sys.argv[1]
outfile = sys.argv[2]
mapfile = sys.argv[3]
ncbi_dir = "ncbi_records/"
newMap(infile, outfile, mapfile, ncbi_dir, amine_receptors)

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	