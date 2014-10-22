# SJS 9/6/14, 10/22/14.
# Run GPCRHMM on sequences, after basic Blast filtering to remove ambiguous, low-quality, and no-corresponding-nuc-record sequences.
# Script ensures that each sequence is a GPCR (with very conservative thresholds) and collects GPCRHMM record if sequence retained.
# This script must be run on Linux to use GPCRHMM locally for GPCR-validation!!

from Bio import AlignIO
import sys
import re
import os
import subprocess
import urllib
import urllib2

GPCR_GLOBAL = 10
GPCR_LOCAL = 10


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
    os.remove("in.fasta")
    os.remove("out.fasta")
    get_scores = re.search('^id\s+(-*\d*\.*\d*)\s+(-*\d*\.*\d*)\s+(\w+)\s*', line)
    assert(get_scores), "Couldn't parse the gpcrhmm output file."
    global_score = get_scores.group(1)
    local_score = get_scores.group(2)
    pred = get_scores.group(3)
    if pred == 'GPCR' and float(global_score) >= GPCR_GLOBAL and float(local_score) >= GPCR_LOCAL:
        return True
    else:
        return False
    

    

def collect_gpcrhmm_structure(protseq, outfile):
    ''' Query the gpcrhmm server (as the released executable can't do this) for the structure/domains. '''

    url1 = 'http://gpcrhmm.sbc.su.se/cgi-bin/predict'
    values = {'protseq':protseq, 'format':'plp'}
    
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
     
        with open(outfile,'w') as openf:
            outf.write(page2)
               

def main():
  
    # Input arguments
    assert(len(sys.argv) == 3), "Usage: python run_process_blast.py <input_protein_file> <input_nucleotide_file> <rdir>. <input_protein_file> is a FASTA sequence file containing filtered protein blast results, and <input_nucleotide_file> is the corresponding FASTA file for nucleotide records. <rdir> is the final output directory."
    naive_protein     = sys.argv[1]  
    naive_nucleotide  = sys.argv[2]
    protein_outfile   = "protein_records_gpcr.fasta"
    nuc_outfile       = "nucleotide_records_gpcr.fasta"
    gpcrhmm_directory = rdir + "/gpcrhmm_records/"
    os.mkdir(gpcrhmm_directory)

    
    # Loop over sequences in <input_protein_file> and validate it as a GPCR. If validated, collect its structure from GPCRHMM server.
    naive_protein = AlignIO.read(naive_protein, "fasta")
    retained_ids = []
    for rec in naive_protein:
        protseq = str(rec.seq)
        fullid  = str(rec.id)
        protid  = fullid.split("_")[0] + "_" + fullid.split("_")[1]
        
        if ensure_gpcr( protseq ):
            
            # This id should be saved.
            retained_ids.append(fullid)
            
            # Retrieve the gpcrhmm record from their server
            collect_gpcrhmm_structure(protseq, gpcrhmm_directory+protid+".txt")
            assert(os.path.exists(gpcrhmm_directory+prot_id+".txt")), "Couldn't get gpcrhmm structural information."
             
            # Save the GPCR-validated protein sequences
            with open(protein_outfile, 'a') as prot_handle:
                prot_handle.write(">" + fullid + "\n" + protseq + "\n")


    # Create new nucleotide sequence fasta file for the GPCR-validated sequences only.
    naive_nuc = AlignIO.read(naive_nucleotide, "fasta")
    nuc_handle = open(nuc_outfile, "w")
    for rec in naive_nuc:
        if str(rec.id) in retained_ids:
            nuc_handle.write(">" + str(rec.id) + "\n" + str(rec.seq) + "\n")
            
        
        
main()