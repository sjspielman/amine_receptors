# SJS. Bootstrap an alignment with fasttree
# USAGE: python bootstrap.py <alignment_file> <datatype> <bootnum> 
### datatype is either protein or dna. gap_exclude is max percentage of gaps allowed in a column to be used. numBoot is number of bootstraps.

from random import randint
import math
import sys
import subprocess
import os
from Bio import AlignIO


class Bootstrap:
	def __init__(self, **kwargs)
			
		self.seqfile = kwargs.get(seqfile, '')		
		assert (os.path.exists(infile)), "Provided input file does not exist."
		self.bootnum = kwargs.get(bootnum, 100)
		self.percent = kwargs.get(percent, 0.90)
		self.datatype = kwargs.get(datatype, '')
		assert(self.datatype=='protein' or self.datatype=='dna'), "\n\nAccepted data types are protein or dna."
		self.cluster = kwargs.get(cluster) # True if on cluster. False if on my computer.
		
		self.seqs = AlignIO.read(infile, 'fasta')
		self.numseq = len(self.seqs)
		self.alnlen = len(self.seqs[0].seq)
		self.limit = math.ceil(percent*self.numseq)
		self.badcol = []
		
		# Bootstrap alignment(s) go to alnfile, bootstrap tree(s) go to treefile
		self.alnfile = 'bootaln.txt'
		self.treefile = 'boottree.txt'

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
		self.cullGap()
		self.makeBootAlignments()
		if self.datatype == 'protein':
			arg = '-wag'
		elif self.datatype == 'dna':
			arg = 'gtr'
		if self.cluster:
			BuildTree='/share/apps/fasttree-2.1.3/FastTree '+arg+' -fastest -nosupport -quiet -n '+str(self.bootnum)+' '+self.alnfile+' > '+self.treefile
		else:
			BuildTree='FastTree '+arg+' -fastest -nosupport -quiet -n '+str(self.bootnum)+' '+self.alnfile+' > '+self.treefile
		runtree=subprocess.call(str(BuildTree), shell='True')	
		assert (runtree == 0), "FastTree fail"

######## INPUT ARGUMENTS ###########
infile = sys.argv[1]
dtype = sys.argv[2]
n = int(sys.argv[3])
cluster_run = int(sys.argv[4]) ## 0 or 1

######## RUN THE BOOTSTRAP #########
boot = Bootstrap(seqfile = infile, datatype = dtype, bootnum = n, cluster = cluster_run )
boot.buildBootTrees()
