# SJS 4/6/14. Basic script to create a map of taxa names, since the names returned have a lot of forbidden characters and tree software freaks out.
#Usage: python unMap.py <infile> <outfile> <mapfile>
## Infile is alignment file with mapped taxa names.
## Outfile is where we will write the alignment with original ids.
## Mapfile is a tdf where column one is original name, column two is map (eg, taxon5) name

from Bio import AlignIO
import sys, re

###########################################################################################################
def readMap(mapfile):
	''' Make a dictionary of mapped names in form {taxon1: histamine_whatever_original_name}. '''
	map_dict = {}
	map = open(mapfile, 'r')
	maplines = map.readlines()
	map.close()
	for line in maplines:
		find = re.search("^(taxon\d+)\t(.+)$", line)
		assert(find), "Bad map search."
		map_dict[find.group(1)] = find.group(2)
	return map_dict
###########################################################################################################

###########################################################################################################
def writeAln(infile, outfile):
	''' Write out alignment with the original names. '''
	records = AlignIO.read(infile, 'fasta')
	outaln = open(outfile, 'w')
	for record in records:
		name = str(record.id)
		try:
			original = map_dict[name]
		except KeyError:
			# presumably this can happen from seed sequences in the alignment (as in pfam gpcrs for UPP)
			print name
			continue
		outaln.write(">"+original+"\n"+str(record.seq)+"\n")
	outaln.close()
###########################################################################################################
###########################################################################################################

infile = sys.argv[1]
outfile = sys.argv[2]
mapfile = sys.argv[3]

map_dict = readMap(mapfile)
writeAln(infile, outfile)














