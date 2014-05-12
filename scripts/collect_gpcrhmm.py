## SJS 5/12/14. 
# Simple script to query gpcrhmm for sequences in a given alignment/sequence file.
# Places all returned files into an output directory named like the sequence file, ish. 


import urllib2
import re
import sys
from Bio import SeqIO

infile = sys.argv[1]

name = infile.split('.')[0]
outdir = 'gpcrhmm_'+name+'/'

records = list(SeqIO.parse(infile, 'fasta'))
for entry in records:
	
	id = str(entry.id)
	seq = str(entry.seq) 
	print id
	
	url1 = 'http://gpcrhmm.sbc.su.se/cgi-bin/predict'
	values = {'protseq':seq, 'format':'plp'}
	
	data = urllib2.urlencode(values)
	req = urllib2.Request(url1, data)
	response1 = urllib2.urlopen(req)
	page1=str(response1.read())
	find_link=re.search('tmp/(\w\w\w\w)', page1)
	if find_link:
		link=str(find_link.group(1))
		url2='http://gpcrhmm.sbc.su.se/tmp/'+link+'.plp'
		response2 = urllib2.urlopen(url2)
		page2=str(response2.read())
	
		outfile=open(outdir+id+'.txt','w')
		outfile.write(page2)
		outfile.close()