##### SJS 2/26/14. Remove duplicates and seqs with too many ambig (>=5%) from HRH blast search.
# Usage: python process_raw_seqs.py <infile> <outfile> <percent of ambiguities we'll accept>


from Bio import SeqIO
import sys


def killDuplicates(records):
	'''Remove duplicate sequences'''
	
	repIDs=0
	repSeqs=0
	count=0
	for ref_record in records:
		w_count = count+1
		while w_count < len(records):
			test_record = records[w_count]
			if ref_record.id != test_record.id:
				# If sequences are the same but different ids, keep only one of the records
				if str(ref_record.seq) == str(test_record.seq):
					repSeqs+=1
					removed_rec = records.pop(w_count)
					w_count-=1
			#If same ID and same seq, remove
			else:
				if str(ref_record.seq) == str(test_record.seq):
					removed_rec = records.pop(w_count)
					w_count-=1
					repIDs+=1
				else:
					print str(ref_record.id)+'\t'+str(test_record.id)+'\n'+str(ref_record.seq)+'\n'
			w_count += 1
		count+=1
	print "duplicate ids:",repIDs
	print "dup seqs with unique ids:",repSeqs

def cullAmbig(records, ambig_percent):
	''' Remove any sequences with >x% ambiguities '''
	bad_amino = ["B", "X", "Z"]
	total_removed = 0
	
	count=0
	for record in records:
		seq = str(record.seq)
		seqlen = float(len(seq))
		badness=0
		for bad in bad_amino:
			badness += seq.count(bad)
		if float(badness)/seqlen > ambig_percent :
			records.pop(count)
			total_removed+=1
		count+=1
	print "seqs with too many ambiguities:", total_removed
################################################################################################
################################################################################################
	
if len(sys.argv) != 4:
	print "Usage: python process_raw_seqs.py <infile> <outfile> <ambigpercent>. (last arg should be decimal float)"
	sys.exit()
	
infile = sys.argv[1]
outfile = sys.argv[2]
ambig_percent = float(sys.argv[3])

records=list(SeqIO.parse(open(infile, 'rU'), "fasta"))
print "original number of seqs:",len(records)
killDuplicates(records)
cullAmbig(records, ambig_percent)
print "final number of seqs:",len(records)

out = open(outfile, 'w')
for record in records:
	out.write(">"+str(record.id)+'\n'+str(record.seq)+'\n')
out.close()