##### Python script used to color certain taxa with figtree. The figtree_colors dictionary was manually created because android color codes are impossible to google, oddly.
import re
import dendropy

figtree_colors = {"lightgreen":"[&!color=#-15483621]", "lightblue":"[&!color=#-14569729]", "darkblue":"[&!color=#-15390531]", "red":"[&!color=#-3400922]", "yellow":"[&!color=#-256]", "orange":"[&!color=#-32768]", "darkgreen":"[&!color=#-13399478]", "cyan":"[&!color=#-15532057]", "rose":"[&!color=#-29506]", "neongreen":"[&!color=#-15008000]", "brown":"[&!color=#-6724045]", "magenta":"[&!color=#-65289]", "purple":"[&!color=#-9363267]"}
receptor_colors = {"histamine": "lightgreen", "trace": "cyan", "hydroxy": "lightblue", "dopamine": "darkblue", "adrenergic": "red", "muscar": "orange"}

### OLD: converted from newick to nexus format for easier coloring ###
infile = "/Users/sjspielman/Desktop/amine_processed_decent2.tre"
outfile = "/Users/sjspielman/Desktop/amine_processed_decent2.nex"
tree = dendropy.Tree.get_from_path(infile, "newick")
tree.write_to_path(outfile, "nexus")

infile = "/Users/sjspielman/Desktop/amine_processed_decent2.nex"
treeh = open(infile, 'rU')
tree = treeh.readlines()
treeh.close()

outfile = "/Users/sjspielman/Desktop/amine_processed_decent2_color.nex"
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
			if entry.upper() in line.upper():
				line = line.rstrip()
				newline = line+figtree_colors[receptor_colors[entry]]+'\n'
				outh.write(newline)
				found=True
		if not found:
			outh.write(line)
	index+=1

for i in range(index, len(tree)):
	outh.write(tree[i])
outh.close()