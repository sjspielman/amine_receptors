##### Python script used to color certain taxa with figtree. The figtree_colors dictionary was manually created because android color codes are impossible to google, oddly.
# USAGE: python colorTree.py infile outfile mapfile
### infile is the newick tree
### outfile is a NEXUS tree file (easier for coloring)
### mapfile is tdf in which column one is original name, column two is map (eg, taxon5) name

import re
import sys
import dendropy

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


figtree_colors = {"lightgreen":"[&!color=#-15483621]", "lightblue":"[&!color=#-14569729]", "darkblue":"[&!color=#-15390531]", "red":"[&!color=#-3400922]", "yellow":"[&!color=#-256]", "orange":"[&!color=#-32768]", "darkgreen":"[&!color=#-13399478]", "cyan":"[&!color=#-15532057]", "rose":"[&!color=#-29506]", "neongreen":"[&!color=#-15008000]", "brown":"[&!color=#-6724045]", "magenta":"[&!color=#-65289]", "purple":"[&!color=#-9363267]"}
receptor_colors = {"histamine": "lightgreen", "trace": "cyan", "hydroxy": "lightblue", "dopamine": "darkblue", "adrenergic": "red", "muscar": "orange"}


infile = sys.argv[1]
outfile = sys.argv[2]
mapfile = sys.argv[3]
map_dict = readMap(mapfile)
tempfile = 'temp.tre' # for non-color nexus tree

# Convert to nexus for easier color parsing
tree = dendropy.Tree.get_from_path(infile, "newick")
tree.write_to_path(tempfile, "nexus")

# Read in nexus tree file
treeh = open(tempfile, 'rU')
tree = treeh.readlines()
treeh.close()

# Write colored tree
outh = open(outfile, 'w')

###### Read until we reach the word "END" for the first time, signifying the end of the taxa block.
found=True
index = 0
for line in tree:
	if "END;" in line:
		break
	else:
		found=False
		for entry in receptor_colors:
			find = re.search('(taxon\d+)', line)
			if find:
				tree_taxon = find.group(1)
				true_taxon = map_dict[tree_taxon]
				print "true", true_taxon.upper()
			
				if entry.upper() in true_taxon.upper():
					line = line.rstrip()
					newline = line+figtree_colors[receptor_colors[entry]]+'\n'
					outh.write(newline)
					found=True
				else:
					print true_taxon				
					
		if not found:
			outh.write(line)
	index+=1

for i in range(index, len(tree)):
	outh.write(tree[i])
outh.close()