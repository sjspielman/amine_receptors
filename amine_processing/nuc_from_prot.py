# SJS 4/16/14.
# Script to retrieve nucleotide ids/seqs from NCBI given all the protein ids.



from Bio import Entrez,SeqIO,Seq
from Bio.Alphabet import generic_dna
Entrez.email = "stephanie.spielman@gmail.com"
import re
import sys
import os
# Regular expression for a CDS ID
cds_match = re.compile("\w\w_\d+\.*\d*")

######################################################################
def matchSeqs(seq1, seq2):
	''' Seqs must be strings of the same length. Just go over and make sure they are identical.'''
	ok = True
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			ok = False
			print i, seq1[i], seq2[i]
			break
	return ok
######################################################################

######################################################################
def fetchRecord(id, type):
	''' Fetch a record from ncbi '''
	assert(type == 'protein' or type == 'nucleotide'), "Invalid ncbi db specification."
	fetch = Entrez.efetch(db=type,  id=id, rettype="gb", retmode="text")
	record = SeqIO.read(fetch, 'gb')
	return record
######################################################################	
		
######################################################################		
def getCDS(prot_ncbi, cds_match):
	'''Retrieve CDS info (id and location) from a given protein record ''' 
	feat = prot_ncbi.features
	cds_id = ''
	for entry in feat:
		if entry.type=='CDS':
			codedby = entry.qualifiers['coded_by'][0]
			
			cds_id = codedby.split(":")[0]
			assert(cds_match.match(cds_id)), "No coding sequence id was found."
			
			cds_location = codedby.split(":")[1]
			
			# Draft genome(s) have some funkiness going on w/ loction.
			if "<" in cds_location or ">" in cds_location:
				cds_location = re.sub(r"<|>", "", cds_location)
				cds_end = int(cds_location.split("..")[1])   # Draft genomes apparently have no stop codons.
			else:
				cds_end = int(cds_location.split("..")[1]) - 3 # -3 to remove the stop codon
			cds_start = int(cds_location.split("..")[0]) - 1 # -1 for indexing
			#print cds_start, cds_end, cds_length, lenprot
			
	return (cds_id, cds_start, cds_end)
######################################################################
	
######################################################################
def main(infile, aa_outfile, nuc_outfile, map_outfile, protid_regexp):
	''' Given a fasta protein sequence file, get corresponding nucleotide data.
		Additionally save all information to a map file for future reference.
	'''
	
	# Open outfiles for writing
	out_aa = open(aa_outfile, 'w')
	out_nuc = open(nuc_outfile, 'w')
	out_map = open(map_outfile, 'w')
	
	seqs = list(SeqIO.parse(infile, 'fasta'))
	for record in seqs:
		# Retrieve the protein id from its fasta sequence record
		full_id = str(record.id)
		find_prot_id = re.search(protid_regexp, full_id)
		assert(find_prot_id), "Couldn't get protid :/"
		prot_id = find_prot_id.group(1)
		protseq = str(record.seq)
		lenprot = len(protseq)
		
		# Retrieve protein NCBI record
		record = fetchRecord(prot_id, "protein")
		
		# Obtain cds info and double check it
		cds_id, cds_start, cds_end = getCDS(record, cds_match)
		cds_length = cds_end - cds_start
		assert(cds_length%3 == 0), "Nucleotide length not a multiple of 3."
		assert(lenprot == cds_length/3), "Nucleotide length /3 is not the same as the protein length."
		
		
		# Grab the CDS record from NCBI and make sure its translation matches original protseq
		nuc_ncbi = fetchRecord(cds_id, "nucleotide")
		nucseq = str(nuc_ncbi.seq)[cds_start:cds_end]	
		nuc2prot = Seq.Seq(nucseq, generic_dna).translate()
		assert ( matchSeqs(nuc2prot, protseq) ), "Translation does not match!"
	

		# Finally, save the final protein and nucleotide files
		final_id = prot_id+'_'+cds_id
		print final_id
		out_aa.write('>'+final_id+'\n'+protseq+'\n')
		out_nuc.write('>'+final_id+'\n'+nucseq+'\n')
		out_map.write(final_id + '\t' + full_id + '\t' + str(record.description) + '\n')
	
	out_aa.close()
	out_nuc.close()
	out_map.close()
######################################################################
		


# Aaaand go.
aa_infile = sys.argv[1]
aa_outfile = sys.argv[2]
nuc_outfile = sys.argv[3]
map_outfile = sys.argv[4]

protid_regexp = "^([XN]P_\d+\.\d)_.+$"
main(aa_infile, aa_outfile, nuc_outfile, map_outfile, protid_regexp)








