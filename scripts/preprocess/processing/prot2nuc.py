# Not yet used, but written to get corresponding nucleotide sequences from the retrieved protein sequences

from Bio import Entrez,SeqIO
import re
######## THESE FILES ARE SO DISGUSTING RANT RANT RANT ###########

Entrez.email = "stephanie.spielman@gmail.com"


cds_match = re.compile("\w+_*\d+\.*\d*")

# Retrieve the coding sequence id given a protein id
prot_id = "ELK28929.1"
prot = Entrez.efetch(db="protein",  id=prot_id, rettype="gb", retmode="text")
record = SeqIO.read(prot, 'gb')
feat = record.features
cds_id = ''
for entry in feat: #loop over the qualifiers
	if entry.type=='CDS':
		cds_id = entry.qualifiers['coded_by'][0]
		cds_id = cds_id.split(":")[0]
		if cds_id.count('(') > 0:
			cds_id = cds_id.split("(")[1]
		
print cds_id
assert(cds_match.match(cds_id)), "No coding sequence id was found."

# Retrieve that coding sequence record
nuc = Entrez.efetch(db="nucleotide",  id=cds_id, rettype="gb", retmode="text")
record = SeqIO.read(nuc,'gb')
feat = record.features
found=False
for entry in feat:
	if entry.type == 'CDS':
		found=True
		# Double check that this nuc record maps back to our original protein id
		if prot_id ==entry.qualifiers['protein_id'][0]:
			print ">"+str(record.id)+'\n'+str(record.seq)+'\n'
		else:
			"poor mapping for prot:",prot_id,"and nuc:",cds_id
if not found:
	print "bad nuc?"