# SJS 9/6/14.
# This script processes the Blast output ids. Retains only sequences for which..
## 1. <1% of the protein sequence is ambiguous, 
## 2. NCBI annotations does not indicate low-quality or similar
## 3. A correct, corresponding nucleotide record in NCBI exists.
## Nucleotide and protein NCBI records for all retained sequences are saved.
## Output files "protein_records_naive.fasta" and "nucleotide_records_naive.fasta" are saved. Note that these files have the same ids!

from Bio import SeqIO, Seq, Entrez
from Bio.Alphabet import generic_dna
import sys
import re
import os
import subprocess
import urllib
import urllib2

Entrez.email = "stephanie.spielman@gmail.com"
BAD_TAGS = ['LOW_QUALITY', 'PSEUDOGENE', 'PARTIAL', 'LOW QUALITY']
RESTRICT_TAXA = 'Vertebrata'
AMBIG_AMINO = ["B", "X", "Z"]
AMBIG_THRESHOLD = 0.01
CDS_MATCH = re.compile("\w\w_\d+\.*\d*")


    
def read_ids(infile):
    ''' Read in the file containing NCBI ids returned from Blast search. Save ids to a set so that we have no duplicates.'''
    id_set = set()
    file = open(infile, 'r')
    lines = file.readlines()
    file.close()
    for line in lines:
        id_set.add(line.strip())
    return list(id_set)
   

def fetch_ncbi_record(id, db, ncbi_dir):
    ''' Fetch an NCBI record if not already fetched and saved. '''
    
    if not os.path.exists(ncbi_dir + id + ".txt"):
        no_response = True
        while no_response:
            try:
                fetch = Entrez.efetch(db=db, rettype="gb", retmode="text", id=id)
                no_response = False
            except:
                pass
        record = SeqIO.read(fetch, "gb")
    else:
        record = SeqIO.read(ncbi_dir + id + ".txt", "gb")
    return record

 
        
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



def keep_ambiguity(sequence):
    ''' Return True if ambiguities make up <=ambig_threshold (this should be a decimal) of the sequence.'''
    total_ambig = 0
    for ambig in AMBIG_AMINO:
        total_ambig += sequence.count(ambig)
    if float(total_ambig)/len(sequence) > AMBIG_THRESHOLD:
        return False
    else:
        return True

    
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

    
def retrieve_nuc_record(prot_record, protseq, ncbi_dir):

    # Obtain CDS info and double check it
    cds_id, cds_start, cds_end = get_CDS_info(prot_record)
    cds_length = cds_end - cds_start
    if cds_length%3 == 0 and len(protseq) == cds_length/3:
        
        # Query NCBI for the nucleotide record
        nuc_record = fetch_ncbi_record(cds_id, "nucleotide", ncbi_dir)
        
        # Check that nucleotide translation is the same as the protein sequence
        nucseq = str(nuc_record.seq)[cds_start:cds_end]    
        protseq_from_nuc = str( Seq.Seq(nucseq, generic_dna).translate() )
        if protseq_from_nuc == protseq:
            return nuc_record, cds_id, nucseq
        else:
            return None, None, None
    else:
        return None, None, None
       


def main():
  
    # Input arguments
    assert(len(sys.argv) == 3), "Usage: python run_process_blast.py <input_file> <rdir>. Input file should contain list of ncbi protein ids to process. rdir is return directory for output gpcrhmm, ncbi files."
    infile = sys.argv[1]  
    rdir = sys.argv[2]
    protein_outfile   = "protein_records_naive.fasta"
    nuc_outfile       = "nucleotide_records_naive.fasta"
    ncbi_directory    = rdir + "/ncbi_records/"
    os.mkdir(ncbi_directory)
        
    # Read in ids, while removing duplicates
    print "Reading in ids and removing duplicates"
    raw_ids = read_ids(infile)
    
    # Loop over ids and decide whether to keep each one or not.
    for prot_id in raw_ids:
        print "Processing", prot_id
        
        # Retreive the NCBI record
        prot_record = fetch_ncbi_record(prot_id, "protein", ncbi_directory)
        protseq = str(prot_record.seq)
        
        # Determine if we should keep based on taxonomy, description, and ambiguities
        print "Assessing taxonomy, description, ambiguities"
        if keep_taxonomy_description(prot_record) and keep_ambiguity(protseq):
                        
            # Fetch the sequence's corresponding nucleotide record from NCBI, and double-check that it plays nicely with the protein record.
            print "Fetching nucleotide record and comparing to protein record"
            nuc_record, nuc_id, nucseq = retrieve_nuc_record(prot_record, protseq, ncbi_directory)
            
            # If fetching the nucleotide record went ok, we can proceed to save both the protein and the nucleotide NCBI records.
            if nuc_record is not None:
                                 
                # Save the NCBI records
                SeqIO.write(prot_record, ncbi_directory + prot_id+".txt", "gb")
                SeqIO.write(nuc_record, ncbi_directory + nuc_id+".txt", "gb")
                
                # Save the sequences
                fullid = prot_id + "_" + nuc_id 
                with open(protein_outfile, 'a') as prot_handle:
                    prot_handle.write(">" + fullid + "\n" + protseq + "\n")
                with open(nuc_outfile, 'a') as nuc_handle:
                    nuc_handle.write(">" + fullid + "\n" + nucseq + "\n")
                
            else:
               print "We are not saving", prot_id, "due to poor nucleotide record."
               
        else:
            print "We are not saving", prot_id, "due to taxonomy/sequence quality/ambiguities."
            
        
main()