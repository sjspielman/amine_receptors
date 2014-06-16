## SJS 5/12/14. 
# Simple script to query gpcrhmm for sequences in a given alignment/sequence file.
# Places all returned files into an output directory named like the sequence file, ish. 

import urllib, urllib2
import re
import os
import sys
import re
from Bio import SeqIO

assert(len(sys.argv)==3), "Usage: python collect_gpcrhmm.py <infile> <gpcr_dir>. infile should be a fasta file of amino acid sequences."

infile = sys.argv[1]
outdir = sys.argv[2]
if not os.path.exists(outdir):
    os.mkdir(outdir)

records = list(SeqIO.parse(infile, 'fasta'))
for entry in records:
    
    id = str(entry.id)
    outname = outdir+id+'.txt'
    if os.path.exists(outname) and os.path.getsize(outname)!=0:
        continue
        
    seq = str(entry.seq)
    seq = seq.replace('-', '') # remove gaps. relevant if file is alignment.


    url1 = 'http://gpcrhmm.sbc.su.se/cgi-bin/predict'
    values = {'protseq':seq, 'format':'plp'}
    
    data = urllib.urlencode(values)
    req = urllib2.Request(url1, data)
    response1 = urllib2.urlopen(req)
    page1=str(response1.read())
    find_link=re.search('tmp/(\w\w\w\w)', page1)
    if find_link:
        link=str(find_link.group(1))
        url2='http://gpcrhmm.sbc.su.se/tmp/'+link+'.plp'
        response2 = urllib2.urlopen(url2)
        page2=str(response2.read())
    
        outfile=open(outname,'w')
        outfile.write(page2)
        outfile.close()