# SJS 10/19/14.
# Script to obtain G-protein binding information. A nice tiny analysis for doing something with the phylogeny.

import urllib, urllib2
from Bio import AlignIO
import re, os

    
def query_predcouple(protid, protseq, outfile):
    ''' Query the PRED_COUPLE2.0 server to obtain G-protein binding information (if not already done and save).
        Save the server result in case.
        Then parse the result to obtain the 4 probabilities of G-protein binding.
        Returns a dictionary of the G-proteins and their probabilities. Note that the PRED_COUPLE people suggest 0.3 as a cutoff (>0.3 = binds).
    '''
    
    # Obtain prediction
    if not os.path.exists(outfile):
        url = 'http://athina.biol.uoa.gr/cgi-bin/bioinformatics/PRED-COUPLE2/pc2++.cgi'   
        data = urllib.urlencode( {'seq':protseq} )
        req = urllib2.Request(url, data)
        no_response = True
        while no_response:
            try:    
                response = urllib2.urlopen(req)
                no_response = False
            except:
                pass        
        result = str( response.read() )
        with open(outfile, 'w') as outf:
            outf.write(result)
    else:
        with open(outfile, 'r') as inf:
            result = inf.read()
    '''
    # Parse prediction
    binding_dict = {}
    regexp = "<i>(G[\w//]+) - </i>\s+<strong><font color=\w+>(\d\.\d+)"
    find_binding = re.findall(regexp, result)
    if len(find_binding) != 4:
        print "bad parsing!!"
        print find_binding
    else:
        for entry in find_binding:
            binding_dict[entry[0]] = entry[1]
            
    return binding_dict
    '''
    
aln = AlignIO.read("/Users/sjspielman/Research/amine_receptors/analysis/alignments/protein_aln.fasta", "fasta")
outdir = "/Users/sjspielman/Research/amine_receptors/analysis/predcouple_records/"
full_dict = {}
for rec in aln:
    seq = str(rec.seq).replace("-","")
    full = str(rec.id).split("_")
    id = full[0] + "_" + full[1]
    print id
    outfile = outdir + id + ".txt"
    query_predcouple(id, seq, outfile)

        




   






