# Script to assess size of each domain segment using gpcrhmm results. 
# I = inner
# O = outer
# M = membrane


import re
import os
import sys
from Bio import SeqIO


amine_subtypes = {'DOPAMINE':             re.compile('.+(D\d|\dD).+'),
			      'HISTAMINE':            re.compile('.+(H\d).+'), 
				  'TRACE':                re.compile('.+RECEPTOR ([1-9]\w*).+'),
				  'ADREN':                re.compile('.+(\w+[- ][1-9]\w*).+'),
				  'ACETYLC':              re.compile('.+RECEPTOR (M\d\w*).+'),
				  'MUSCARINIC':           re.compile('.+MUSCARINIC (\d\w*).+'),
				  'MUSCARINIC RECEPTOR':  re.compile('.+ MUSCARINIC RECEPTOR (\d\w*) SUBTYPE.+'),
				  'SEROTONIN 5-HT':       re.compile('.+SEROTONIN 5-HT(\d).+'),
                  'SEROTONIN':            re.compile('.+RECEPTOR ([1-9]\w*).+'),
				  '5-HT':                 re.compile('.+5-HT(\d).+'), 
				  '5-HYDROXYTRYPTAMINE':  re.compile('.+RECEPTOR ([1-9]\w*).+')
			     }
amine_wonky = ["UNCHARACTERIZED", "HYPOTHETICAL", "OCTOPAMINE", "PROBABLE", "SPERMIDINE"]

   
def findSubType(prot_id, ncbi_path, dict, wonky):
	''' From a protein id, get some info from ncbi record about subtype.'''
	
	file = ncbi_path + prot_id + '.gb'
	ncbi_record = SeqIO.read(file, 'gb')
	desc = str(ncbi_record.description)
	#for regexps
	desc = 'wompwomp   ' + desc.upper()
	desc = re.sub('\(', '', desc)
	desc = re.sub('\)', '', desc)
	
	subtype = None
	if "UNCHARACTERIZED" in desc or "HYPOTHETICAL" in desc or "OCTOPAMINE" in desc or "PROBABLE" in desc or "SPERMIDINE" in desc:
		subtype = 'unknown'
		return subtype

	for type in dict:
		if type in desc:
			try:
				subtype = dict[type].match(desc).group(1)
				break
			except:
				pass
	if subtype is not None:
		return subtype
	else:
		for wonk in wonky:
			if wonk in desc:
				subtype = 'unknown'
				return subtype
	if subtype is None:
		print "boooooo"
		sys.exit()




def whichDomain(inner, outer, membr, cutoff):
	''' Given probabilities for each domain, which are we likely in? Lenient ambiguity, please! '''
	temp = [ inner, outer, membr ]
	which = temp.index(max(temp))
	if which == 0:
		return 'I'
	elif which == 1:
		return 'O'
	elif which == 2:
		return 'M'
	else:
		raise AssertionError("How in the hell did we get here?!!???!!!!")


##########################################################################################	
def parseGPCRHMM(file, cutoff):
	''' Parse a gpcrhmm file and return a list of domain letters to replace amino acid sequence with. ''' 
	hmm_file = open(file, 'rU')
	hmm_lines = hmm_file.readlines()
	hmm_file.close()
	
	sizeList = []
	previousDomain = 'O' # Start with O since all begin as outer/extracelluar Nterm.
	size = 0
	for line in hmm_lines:
		parseLine=re.search('^(\d+)\s\w\t(.+)\t(.+)\t(.+)\t(.+)$', line)
		if parseLine:
			total = parseLine.group(1)
			inner = float(parseLine.group(2))
			outer =  float(parseLine.group(3)) + float(parseLine.group(5))
			membr = float(parseLine.group(4))
			currentDomain = whichDomain(inner, outer, membr, cutoff)
			if currentDomain == previousDomain:
				size += 1
			else:
				sizeList.append(str(size))
				previousDomain = currentDomain
				size = 1
	sizeList.append(str(size))
	return (total, sizeList)

def readMap(mapfile):
	''' Parse map to dictionary of form {name: ncbi_id}. '''
	map_dict = {}
	map = open(mapfile, 'r')
	maplines = map.readlines()
	map.close()
	for line in maplines:
		find = re.search("^(.+)\t(.+)$", line)
		assert(find), "Bad map search."
		map_dict[find.group(1)] = find.group(2)
	return map_dict	
		
		
def main(infile, map, outfile, ncbi_path, hmm_dir, cutoff):
	''' Get size of each domain using gpcrhmm output.
		infile  = sequence infile
		map     = duh?
		outfile = duh?
		hmm_dir = directory to the gpcrhmm files
		cutoff  = posterior probability threshold for calling domain.
	'''
	aln = list(SeqIO.parse(infile, 'fasta'))
	outf = open(outfile, 'w')
	outf.write('name	type	full	Nterm	M1	ICL1	M2	ECL1	M3	ICL2	M4	ECL2	M5	ICL3	M6	ECL3	M7	Cterm\n')
	for record in aln:
		id = str(record.id)
		id_generic = id.split("_")[0]
		print id
		ncbi = map[id]
		type = findSubType(ncbi, ncbi_path, amine_subtypes, amine_wonky)
		seq = str(record.seq)
		hmmFile = hmm_dir + id + '.txt'
		total, sizeList = parseGPCRHMM(hmmFile, cutoff)
		outf.write(id + '\t' + id_generic + '\t' + type + '\t' + total + '\t' + '\t'.join(sizeList)+'\n')
	outf.close()
##########################################################################################	

path = os.getcwd() +'/'

if len(sys.argv) != 6:
	print "Usage: python getDomainSizes.py <infile> <mapfile> <outfile> <hmm_dir> <cutoff> Paths not needed, but assumes infile and hmm_dir are in same directory."
	sys.exit()

infile = path + sys.argv[1]
mapfile = path + sys.argv[2]
outfile = path + sys.argv[3]
hmm_dir = path + sys.argv[4]
if hmm_dir[-1] != '/':
	hmm_dir += '/'
cutoff = float(sys.argv[5])
ncbi_path = path + 'ncbi_records/'
map = readMap(mapfile)
main(infile, map, outfile, ncbi_path, hmm_dir, cutoff)
















