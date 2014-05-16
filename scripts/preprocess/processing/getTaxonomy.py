# misc crap for entrez

from Bio import Entrez,SeqIO
import re

Entrez.email = "stephanie.spielman@gmail.com"


# Contains all sequences used in phylogeny. This isn't an alignment, though.
seqfile = '/Users/sjspielman/Dropbox/HRH/amine_3_2_14/human_seeds/amine_processed.fasta'


keep_taxa = ['HomoSapiens', 'MusMusculus', 'GallusGallus', 'DanioRerio', 'XenopusLaevis']

outfile = '/Users/sjspielman/Dropbox/HRH/amine_3_2_14/human_seeds/amine_processed_vert.fasta'
out = open(outfile, 'a')

count=0
#removed=252
seqs = SeqIO.parse(seqfile, 'fasta')
for record in seqs:
	#if count < 3463:
	#	count+=1
	#	continue
	print count
	find = re.search("^(\w\w_\d+\.\d)_.+", str(record.id))
	if find:
		prot_id = find.group(1)
		print prot_id
		fetch = Entrez.efetch(db="protein",  id=prot_id, rettype="gb", retmode="text")
		protSeq = SeqIO.read(fetch, 'gb')
		print protSeq.description.upper()
		print protSeq.annotations['comment']
		print
		#if "PREDICTED" in protSeq.description:
		#	print "wompwomp"
		taxonomy =  ("").join(protSeq.annotations["taxonomy"])
		print taxonomy
		if "Vertebrata" in taxonomy:
			print "hurray"
		#	out.write(">"+str(record.id)+"\n"+str(record.seq)+"\n")			
		#else:
		#	removed+=1
		#	print "\t"+str(removed)
	else:
		assert(1==0), "poop"
	count+=1
out.close()









'''
dontwant = ["LOW_QUALITY", "putative", "protein", "related", "precursor", "pseudogene", "like", "predicted", "similar", "uncharacterized", "isoform", "-P", "GF"]

		
seqs = SeqIO.parse(seqfile, 'fasta')
for record in seqs:
	ok = True
	find = re.search("(\w\w_\d+\.\d)_(.+)", str(record.id))
	if find:
		for entry in dontwant:
			if entry in find.group(2):
				ok = False
				break
		if ok:
			# If the annotation isn't gross, grab the taxonomy.
			prot_id = find.group(1)
			fetch = Entrez.efetch(db="protein",  id=prot_id, rettype="gb", retmode="text")
			protSeq = SeqIO.read(fetch, 'gb')
			print protSeq.annotations["taxonomy"][-3:]

'''
'''


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
'''