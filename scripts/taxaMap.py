# SJS 3/13/14. Basic script to create a map of taxa names, since the names returned have a lot of forbidden characters and tree software freaks out.
#Usage: python taxaMap.py <infile> <outfile> <mapfile>
## Infile a fasta sequence file (aligned or unaligned, doesn't matter) with original taxa names.
## Outfile is where we will write the sequences with their new ids
## Mapfile is a tdf where column one is old name, column two is new name

from Bio import SeqIO
import sys


infile = sys.argv[1]
outfile = sys.argv[2]
mapfile = sys.argv[3]

records = list(SeqIO.parse(infile, 'fasta'))

outaln = open(outfile, 'w')
map = open(mapfile, 'w')

counter = 0
for record in records:
		original = str(record.id)
		new = "taxon"+str(counter)
	
		outaln.write(">"+new+"\n"+str(record.seq)+"\n")
		map.write(new+'\t'+original+'\n')
		
		counter+=1

outaln.close()
map.close()