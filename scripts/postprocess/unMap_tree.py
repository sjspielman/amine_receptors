# SJS 4/6/14. Basic script to create a map of taxa names, since the names returned have a lot of forbidden characters and tree software freaks out.
#Usage: python unMap.py <infile> <outfile> <mapfile>
## Infile is alignment file with mapped taxa names.
## Outfile is where we will write the alignment with original ids.
## Mapfile is a tdf where column one is original name, column two is map (eg, taxon5) name

from Bio import AlignIO
import sys
import re
import dendropy

###########################################################################################################
def swapTaxon(tree_string, map_dict, index):
	find = re.search('^(taxon\d+)[,:\)]', tree_string)
	assert(find), "Can't parse taxon name."
	map_name = find.group(1)
	print map_name
	#assert 1==0
	try:
		orig_name = map_dict[map_name]
	except KeyError:
		print "womp womp!!!"
		print map_name
		sys.exit()		
	
	return(orig_name, index + len(map_name))
###########################################################################################################
def parseTree(infile, outfile, map_dict, index = 0):
	inf = open(infile, 'r')
	tree_string = inf.read()
	inf.close()
	tree_string = re.sub(r"\s", "", tree_string)
	treelen = len(tree_string)
	out_tree = ''
	
	while index < treelen:
		print index  #, out_tree
		if tree_string[index] == 't':
			add_string, index = swapTaxon(tree_string[index:], map_dict, index)
			out_tree += add_string
		else:
			out_tree += tree_string[index]
			index += 1	
	
	outf = open(outfile, 'w')
	outf.write(out_tree)
	outf.close()
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

infile = sys.argv[1]
outfile = sys.argv[2]
mapfile = sys.argv[3]

if len(sys.argv) != 4:
	print "Usage: python unMap.py <infile> <outfile> <mapfile>"
	sys.exit()

map_dict = readMap(mapfile)
parseTree(infile, outfile, map_dict)













