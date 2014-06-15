# SJS 5/16/14.
# Simply fetch NCBI records and save to avoid future fetching. Can simultaneously cull, as desired.

import sys
import re
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "stephanie.spielman@gmail.com"

###########################################################################################################
def cullFromDesc(string):
    ''' Find if has something bad.... excellent documentation, steph. really good.'''
    bad = ['LOW_QUALITY_PROTEIN', 'PSEUDOGENE', 'PARTIAL', 'LOW QUALITY PROTEIN']
    string = string.upper()
    for entry in bad:
        if entry in string:
            return False
    return True
###########################################################################################################    
def getRecordAndCull(prot_id, raw_id, ncbi_dir, restrictClade = False):
    ''' Retrieve NCBI record. 
        We also have the option to save only members of a clade (restrictClade).
        Also also wik, don't save anything crappy. 
    '''

    fetch = Entrez.efetch(db="protein",  id=prot_id, rettype="gb", retmode="text")
    ncbi_record = SeqIO.read(fetch, 'gb')
    keep = cullFromDesc(str(ncbi_record.description))
    if keep:
        if restrictClade:
            fullTax = ncbi_record.annotations["taxonomy"]
            if restrictClade not in fullTax:
                print "DISCARDING ", raw_id," - NOT A "+restrictClade.upper()
                return False
        outfile = ncbi_dir + prot_id + '.gb'
        SeqIO.write(ncbi_record, outfile, 'gb')    
    else:
        print "DISCARDING ", raw_id," - CRAPPY PROTEIN."
        return False
    return True    
###########################################################################################################
def parseInput(args):
    assert(3 <=len(args) <= 5), "Usage: python fetchNCBI.py <infile> <outfile> <ncbi_directory> <clade restriction> . Last 2 args optional."
    infile = args[1]
    outfile = args[2]
    try:
        ncbi_dir = args[3]
    except:
        ncbi_dir = os.getcwd() + '/ncbi_records/'
    if not os.path.exists(ncbi_dir):
        os.mkdir(ncbi_dir)    
    try:    
        restrictClade = args[4]
    except:
        restrictClade = False
    return(infile, outfile, ncbi_dir, restrictClade)
###########################################################################################################
def main(args):
    ''' Read in records from a sequence file. Retrieve its ncbi record if not yet retrieved.
        Will exclude low quality and seqs which fall outside of a restricted clade, if specified.
        Save fasta file for sequences which proceed in pipeline.
    '''
    
    # Parse input
    infile, outfile, ncbi_dir, restrictClade = parseInput(sys.argv)
    
    # read in records. for each record, fetch the record if it doesn't exist yet. records are named according to ids.
    records = list(SeqIO.parse(infile, 'fasta'))
    outhandle = open(outfile, 'w')
    for record in records:
        raw_id = str(record.id)
        seq = str(record.seq)
        find = re.search("(\wP_\d+\.\d)_(.+)$", raw_id)
        assert(find), "Couldn't get id search."
        prot_id = find.group(1)

        # if don't have already, go get it.
        if not (os.path.exists(ncbi_dir + prot_id + '.gb')):
            keep = getRecordAndCull(prot_id, raw_id, ncbi_dir, restrictClade)
            if keep:
                outhandle.write('>'+raw_id+'\n'+seq+'\n')
    outhandle.close()    
###########################################################################################################


main(sys.argv)







