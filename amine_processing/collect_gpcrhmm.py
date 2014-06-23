## SJS 6/23/14.


import urllib, urllib2
import re
import os
import sys
import subprocess
import shutil
from Bio import SeqIO

assert(len(sys.argv)==2), "Usage: python collect_gpcrhmm.py <infile> . infile should be a fasta file of amino acid sequences."

infile = sys.argv[1]
path = "/".join( infile.split('/')[:-1] )
name_raw = infile.split('/')[-1]

name = name_raw.split('.')[0]
olddir = '/Users/sjspielman/Dropbox/Amine/AmineNew/amine_gpcrhmm/'
newdir = '/Users/sjspielman/Dropbox/Amine/AmineNew/amine_gpcrhmm_new/'


records = list(SeqIO.parse(infile, 'fasta'))
for i in range(len(records)):
    
    print i
    entry = records[i]
    
    id = str(entry.id)
    seq = str(entry.seq)
    find = re.search('^([XN]P_\d+\.\d)_', id)
    assert(find), "boooo no find"
    protid = find.group(1)
    existing_file = olddir + protid + '.txt'
    new_file = newdir + id + '.txt'
    
    if os.path.exists(existing_file):
        shutil.copy(existing_file, new_file)
        continue
    else:    
        print "nope", id
        continue
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
    
            outfile=open(new_file,'w')
            outfile.write(page2)
            outfile.close()