# SJS 3/12/14. Bootstrap an alignment with fasttree
# USAGE: python bootstrap.py <alignment_file> <datatype> <gap_exclude> <numBoot> <bootaln_file> <boottree_file>
### datatype is either protein or dna. gap_exclude is max percentage of gaps allowed in a column to be used. numBoot is number of bootstraps.

from random import randint
import math
import sys
import subprocess
from Bio import AlignIO


class Bootstrap():
	def __init__(self, bootnum, percent, infile, alnfile, treefile, datatype):
		self.seqs = AlignIO.read(infile, 'fasta')
		self.numseq = len(self.seqs)
		self.alnlen = len(self.seqs[0].seq)
		self.limit = math.ceil(percent*self.numseq)
		self.badcol = []
		self.bootnum = bootnum
		self.alnfile = alnfile
		self.treefile = treefile
		self.datatype = datatype
		assert(self.datatype=='protein' or self.datatype=='dna'), "\n\nAccepted data types are protein or dna."
	
	def cullGap(self):
		for i in range(self.alnlen):
			findgaps = str(self.seqs[:,i])
			if findgaps.count('-') > self.limit:
				#print findgaps.count('-'), self.limit
				self.badcol.append(i)
		#print len(self.badcol), self.alnlen

	def makeBootAlignments(self):
		outhandle=open(self.alnfile, 'w')
		for i in range(self.bootnum):
			print "bootaln",i
			outhandle.write(' '+str(self.numseq)+' '+str(self.alnlen)+'\n')
			indices = []
			for a in range(self.alnlen):
				colnum = randint(0,self.alnlen-1)
				while colnum in self.badcol:
					colnum = randint(0,self.alnlen-1)
				indices.append(colnum)	
			for s in range(self.numseq):
				newseq=''
				id=self.seqs[s].id
				for a in indices:
					newseq=newseq+self.seqs[s][a]
				outhandle.write(str(id)+'        '+newseq+'\n')
		outhandle.close()

	def buildBootTrees(self):
		#self.cullGap()
		self.makeBootAlignments()
		if self.datatype == 'protein':
			arg = '-wag'
		elif self.datatype == 'dna':
			arg = 'gtr'
		BuildTree='FastTree '+arg+' -fastest -nosupport -quiet -n '+str(self.bootnum)+' '+self.alnfile+' > '+self.treefile
		runtree=subprocess.call(str(BuildTree), shell='True')	

def parse_args():
    parser = argparse.ArgumentParser(prefix_chars='+-', usage='--protein_file <User File>')
    parser.add_argument("-infile", help="A file containing unaligned sequences in FASTA format", required=False, dest="infile", type=str)
    parser.add_argument("-n",dest="threads", type=int, help="Number of processes to use")
    parser.add_argument("-form",dest="form", type=str, help="The infile format (usually FASTA)", default="FASTA")
    parser.add_argument("-bootstraps", help="The number of bootstraps to perform", required=False,
            dest="bootstraps")
    parser.add_argument("-alphabet", help="Whether AAs or NTs are used (protein or nucleotide)", type=str,
            default="protein", required=False, dest="alphabet") ##AA or NT, default is AA
    ## Gap penalization is now hard-coded by default in accordance with the original runs

    return parser.parse_args()


infile = sys.argv[1]
datatype = sys.argv[2]
gap_limit = float(sys.argv[3])
bootaln = sys.argv[4]
boottree = sys.argv[5]


boot = Bootstrap(100, 0.90, infile, bootaln, boottree)
boot.buildBootTrees()
