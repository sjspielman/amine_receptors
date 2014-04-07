# For sequences returned by our blast search, retain only certain ones based on various genbank annotations.

from Bio import Entrez,SeqIO
import re

Entrez.email = "stephanie.spielman@gmail.com"

exclude_desc = ['PREDICTED', 'MODEL', 'WGS', 'INFERRED', 'LOW_QUALITY', 'PUTATIVE', 'RELATED', 'PRECURSOR', 'PSEUDOGENE', 'LIKE', 'UNCHARACTERIZED']



# Blast sequences culled of duplicates and sequences with excessive ambiguity (>5% i think)
home = '/home/sjs3495/AmineReceptors/human_seeds/'


seqfile = 'amine_processed.fasta'
outfile_decent = 'amine_processed_decent.fasta'
outfile_great  = 'amine_processed_great.fasta'


out_decent = open(outfile_decent, 'w')
out_great = open(outfile_great, 'w')

count=0
seqs = SeqIO.parse(seqfile, 'fasta')
for record in seqs:
	good=True
	#if count < 3463:
	#	count+=1
	#	continue
	print count
	find = re.search("^(\w\w_\d+\.\d)_.+", str(record.id))
	if find:
		prot_id = find.group(1)
		fetch = Entrez.efetch(db="protein",  id=prot_id, rettype="gb", retmode="text")
		protSeq = SeqIO.read(fetch, 'gb')
		desc = protSeq.description
		for entry in exclude_desc:
			if entry in desc:
				good=False
				break
		if good:
			desc = re.sub(" ","_",desc)
			desc = re.sub('\[','',desc)
			desc = re.sub('\]','',desc)
			newID = prot_id+"_"+desc
			print newID
			out_decent.write(">"+newID+'\n'+str(record.seq)+'\n')
			if "VALIDATED" in desc or "REVIEWED" in desc:
				print "Great!!!!"
				out_great.write(">"+newID+'\n'+str(record.seq)+'\n')
	else:
		assert(1==0), "poop"
	count+=1
out_decent.close()
out_great.close()