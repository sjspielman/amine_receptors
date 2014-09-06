# SJS 9/6/14.

from Bio import SeqIO, Seq, Entrez
from Bio.Alphabet import generic_dna
import sys
import re
import os
import subprocess
import urllib
import urllib2


BAD_TAGS = ['LOW_QUALITY', 'PSEUDOGENE', 'PARTIAL', 'LOW QUALITY']
RESTRICT_TAXA = 'Vertebrata'
AMBIG_AMINO = ["B", "X", "Z"]
GPCR_GLOBAL = 10
GPCR_LOCAL = 10
CDS_MATCH = re.compile("\w\w_\d+\.*\d*")

protein_outfile   = "protein_records.fasta"
nuc_outfile       = "nucleotide_records.fasta"
NCBI_directory    = "ncbi_records/"
GPCRHMM_directory = "gpcrhmm_records/"


	
def read_ids(infile):
    ''' Read in the file containing NCBI ids returned from Blast search. Save ids to a set so that we have no duplicates.'''
    id_set = set()
    file = open(infile, 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        id_set.add(line.strip())
    return list(id_set)
   
   
        
def keep_taxonomy_description(record):
    ''' Return True if taxonomy and description are acceptable, otherwise False. '''
    taxonomy = (" ").join(record.annotations["taxonomy"])
    if RESTRICT_TAXA in taxonomy:
        desc = str(record.description).upper()
        for tag in BAD_TAGS:
            if tag in desc:
                return False
        return True
    else:
        return False



def keep_ambiguity(sequence, ambig_threshold):
    ''' Return True if ambiguities make up <=ambig_threshold (this should be a decimal) of the sequence.'''
    total_ambig = 0
    for ambig in AMBIG_AMINO:
        total_ambig += sequence.count(ambig)
	if float(total_ambig)/seqlen > ambig_percent :
        return False
    else:
        return True



def ensure_gpcr(sequence):
    ''' Run gpcrhmm to determine if we really have a GPCR here. Use conservative thresholds. ''' 
    
    # Create input fasta file for gpcrhmm, and execute.
    temp = open("in.fasta", "w")
    temp.write(">id\n"+sequence+"\n")
    temp.close()
	run = subprocess.call("perl gpcrhmm/gpcrhmm.pl in.fasta > out.fasta", shell=True)
	assert(run == 0), "gpcrhmm did not properly run and returned with exit code 0."
	
	# Parse the output file
	outf = open("out.fasta", "rU")
	line = outf.readlines()[1]
	outf.close()
	get_scores = re.search('^'+seqid+'\s+(-*\d*\.*\d*)\s+(-*\d*\.*\d*)\s+(\w+)\s*', line)
	assert(get_scores), "Couldn't parse the gpcrhmm output file."
    global_score = get_scores.group(1)
	local_score = get_scores.group(2)
	pred = get_scores.group(3)
	if pred == 'GPCR' and float(global_score) >= GPCR_GLOBAL and float(local_score) >= GPCR_LOCAL:
	    return True
	else:
	    return False

	
def get_CDS_info(prot_record):
    '''Retrieve CDS info (id and location) from a given protein record ''' 
    
    feat = prot_record.features
    cds_id = ''
    for entry in feat:
        if entry.type=='CDS':
            codedby = entry.qualifiers['coded_by'][0]
            cds_id = codedby.split(":")[0]
            assert(CDS_MATCH.match(cds_id)), "No coding sequence id was found."
                
            cds_location = codedby.split(":")[1]
            if "<" in cds_location or ">" in cds_location:       # Draft genome(s) have some bizarre characters and also no stop codons, so have to do some workarounds here.
                cds_location = re.sub(r"<|>", "", cds_location)
                cds_end = int(cds_location.split("..")[1])
            else:
                cds_end = int(cds_location.split("..")[1]) - 3 # -3 to remove the stop codon
            cds_start = int(cds_location.split("..")[0]) - 1 # -1 for indexing
      
    return cds_id, cds_start, cds_end

    
def retrieve_nuc_record(prot_record, protseq):

    # Obtain CDS info and double check it
    cds_id, cds_start, cds_end = get_CDS_info(prot_record)
    if (cds_end - cds_start)%3 == 0 and len(protseq) == cds_length/3):
        
        # Query NCBI for the nucleotide record
        fetch = Entrez.efetch(db="nucleotide", id=cds_id, rettype="gb", retmode="text")
        nuc_record = SeqIO.read(fetch, 'gb')
        
        # Check that nucleotide translation is the same as the protein sequence
        nucseq = str(nuc_record.seq)[cds_start:cds_end]    
        protseq_from_nuc = str( Seq.Seq(nucseq, generic_dna).translate() )
        if protseq_from_nuc == proseq:
            return nuc_record, nucseq
        else:
            return None, None
    else:
        return None, None
       

def collect_gpcrhmm_structure(protseq, outfile):
    ''' Query the gpcrhmm server for the structure/domains. '''

    url1 = 'http://gpcrhmm.sbc.su.se/cgi-bin/predict'
    values = {'protseq':proseq, 'format':'plp'}
    
    data = urllib.urlencode(values)
    req = urllib2.Request(url1, data)
    no_response1 = True
    while no_response1:
        try:    
            response1 = urllib2.urlopen(req)
            no_response1 = False
        except:
            pass        
    page1=str(response1.read())
    find_link=re.search('tmp/(\w\w\w\w)', page1)
    if find_link:
        link=str(find_link.group(1))
        url2='http://gpcrhmm.sbc.su.se/tmp/'+link+'.plp'
        no_response2 = True
        while no_response2:
        try:    
            response2 = urllib2.urlopen(url2)
            no_response2 = False
        except:
            pass     
        page2=str(response2.read())
     
        outf=open(outfile,'w')
        outf.write(page2)
        outf.close()  

               

def main():
  
    # Input argument
    assert(len(sys.argv) == 2), "Usage: python run_process_blast.py <input_file>. Input file should contain list of ncbi protein ids to process."
    infile = sys.argv[1]  
        
    # Read in ids, while removing duplicates
    raw_ids = read_ids(infile)
    
    # Loop over ids and decide whether to keep each one or not.
    for prot_id in raw_ids:
        
        # Retreive the NCBI record
        fetch = Entrez.efetch(db="protein", rettype="gb", retmode="text", id=prot_id)
        prot_record = SeqIO.read(fetch, "gb")
        protseq = str(record.seq)
        
        # Determine if we should keep based on taxonomy, description, and ambiguities
        if keep_taxonomy_description(prot_record) and keep_ambiguity(protseq, ambig_threshold):
            
            # Ensure that the sequence is really a GPCR. Note that we don't include this in the previous line in order to avoid excessive calls to gpcrhmm
            if ensure_gpcr(protseq):
                
                # Fetch the sequence's corresponding nucleotide record from NCBI, and double-check that it plays nicely with the protein record.
                nuc_record, nucseq = retrieve_nuc_record(prot_id, protseq)
                
                # If fetching the nucleotide record went ok, we can proceed to save both the protein and the nucleotide NCBI records.
                if nuc_record is not None:
                    
                    # Retrieve the gpcrhmm record from their server
                    collect_gpcrhmm_structure(protseq, GPCRHMM_directory+prot_id+".txt")
                     
                    # Save the NCBI records
                    SeqIO.write(prot_record, prot_id+".txt", "gb")
                    SeqIO.write(nuc_record, nuc_id+".txt", "gb")
                    
                    # Save the sequences
                    prot_handle = open(protein_outfile, 'a')
                    prot_handle.write(">"+prot_id+"\n"+protseq+"\n")
                    prot_handle.close()
                    nuc_handle = open(nucleotide_outfile, 'a')
                    nuc_handle.write(">"+nuc_id+"\n"+nucseq+"\n")
                    nuc_handle.close()
    
main()