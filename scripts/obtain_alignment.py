#! /usr/bin/env python

##############################################################################
## Script to obtain structurally-curated alignment of amine receptors.
## Provide two corresponding FASTA sequence files (protein, nucleotide) of GPCRs as well as their gpcrhmm records. 
## IMPORTANT: sequence files must use the same IDs for corresponding sequences. NO checking is done for this!!
## Yields the following files:
## 1-2. <protein/nucleotide>_aln_naive.fasta  : Alignments simply created with mafft,       
## 3-4. <protein/nucleotide>_aln.fasta        : Alignments where domains properly align,    
## 5-6. <protein/nucleotide>_aln_masked.fasta : Alignments where domains properly align, but structurally ambiguous columns and residues which do not match the consensus column structure are masked, 
## 7.  domain_alignment.fasta                 : "Alignment" where residues are replaced with domain tags O/I/M/A (outer, inner, membrane, ambiguous)
## 8.  domain_consensus.txt                   : Simple string of the consensus domains 
##
##
##  Written by Stephanie J. Spielman (stephanie.spielman@gmail.com).
##############################################################################



import sys
import re
import os
import subprocess
import shutil
from string import maketrans
import numpy as np
from Bio import AlignIO, SeqIO



################### GLOBALS ####################
GPCR_DOMAIN_DICT   = {'O':0, 'I':1, 'M':2, 'A':3, '-':4}
TRANSLATOR         =  maketrans("".join(map(str,GPCR_DOMAIN_DICT.values())), "".join(GPCR_DOMAIN_DICT.keys())) # to translate integers to domains
GPCR_STRUCTURE     = [0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1, 3] # we use the 3 as a "hack" to terminate the protein.
PERCENT_MISALIGNED = 0.05  # We discard *sequences* where >=threshold of NONGAP positions do not match consensus structure
GPCRHMM_THRESHOLD  = 0.5   # We assign domains if posterior probability is >= threshold. Otherwise, we call it ambiguous.
################################################




def call_mafft(infile, outfile, options=''):
    ''' Infile contains unaligned protein fasta sequences to be sent to outfile. Options can be given to mafft as well.'''
    mafft = "mafft " + options + " " + infile + " > " + outfile
    runmafft = subprocess.call(mafft, shell=True)
    assert(runmafft == 0), "\n\nMafft failed to run."



def seq_to_gpcrhmm(seqfile, gpcrhmm_directory):
    ''' Convert sequences to their corresponding domains. Note that these sequences are not alignments, so we should use python lists.'''
    seqs = list(SeqIO.parse(seqfile, 'fasta'))
    seq_domains = []
    for record in seqs:
        id = str(record.id)
        gpcrhmm_file = gpcrhmm_directory +  id.split("_")[0]+"_" + id.split("_")[1] + '.txt'
        seq_domains.append( parse_gpcrhmm(gpcrhmm_file, len(str(record.seq))) )
    return seq_domains


def parse_gpcrhmm(file, length):
    ''' Parse a gpcrhmm file and return a list of domain letters to replace amino acid sequence with. ''' 
    
    domain_order = 'IOM'
    putative_domains = np.zeros([3, length]) # first row inner, second row outer, third row membrane.
    with open(file, 'rU') as hmm_file:
        i = 0
        for line in hmm_file:
            parse_line=re.search('^\d+\s\w\t(.+)\t(.+)\t(.+)\t(.+)$', line)
            if parse_line:
                putative_domains[0][i] = float(parse_line.group(1))
                putative_domains[1][i] = float(parse_line.group(2)) + float(parse_line.group(4))
                putative_domains[2][i] = float(parse_line.group(3))
                i += 1   
    domains = ''
    for i in range(length):
        temp = putative_domains[:,i]
        if np.max(temp) >= GPCRHMM_THRESHOLD:        
            domains += domain_order[ np.argmax(temp) ]
        else:
            domains += 'A'
    assert(len(domains) == length), "\n\nCould not identify a domain for all residues."
    return domains
    
    
    
    
def aln_to_gpcrhmm(aln, numseq, numcol, domains):
    ''' Convert a protein alignment file to a "structural" alignment. Replace amino acid residues with domain (O, I, M, A).'''
 
    struc_aln = np.zeros( [numseq, numcol], dtype = 'int8')
    num_aa = np.zeros(numseq, dtype = 'int16')    # Save the number of amino acids per sequence (row in struc_aln)
    i=0
    for record in aln:
        id = str(record.id)
        seq = str(record.seq)
        num_aa[i] = len(seq) - seq.count('-')
        domain_sequence = replace_residues(seq, domains[i])
        struc_aln[i] = domain_sequence
        i += 1
    return struc_aln, num_aa



def replace_residues(seq, domains):
    ''' Replace the amino acid string (which includes gaps!!) with the domain letters. '''
    nogap_count = 0
    total_count = 0
    newseq = np.zeros(len(seq), dtype = 'int8')
    newseq[newseq == 0] = 5 
    
    for position in seq:
        if position == '-':
            newseq[total_count] = 4
        else:
            newseq[total_count] = GPCR_DOMAIN_DICT[ domains[nogap_count] ]
            nogap_count += 1
        total_count += 1
    assert(newseq.all() != 5), "\n\nCould not convert sequence to domains."
    return newseq




def build_keep_list(struc_aln, numseq, num_aa):
    
    domain_index = 0
    consensus_domain = 0
    keep_list = np.zeros(numseq) # Each element is a row in the alignment. The resulting numbers represent how many columns in that sequence violate consensus.
    for column in struc_aln.T:
            
        # Determine the column consensus
        domain_index, consensus_domain = determine_consensus(numseq, column, domain_index, consensus_domain)

        # Add 1 to keep_list for those rows which do not match consensus
        masks = [column != 3, column != 4, column != consensus_domain]
        keep_list += reduce(np.logical_and, masks)
    
    # Return boolean array for whether we should keep each sequence.
    keep_list /= num_aa
    return keep_list < PERCENT_MISALIGNED



def determine_consensus(numseq, column, domain_index, consensus_domain):

    # Determine the column consensus
    unique, positions = np.unique(column, return_inverse=True)
    max = unique[np.bincount(positions).argmax()]
    
    # If the most frequent domain in the column occurs for >=50% of sequences (rows), we consider it the consensus domain. Otherwise, we assume that the domain has not changed from the previous column and we do not reassign.
    percent_max = float(np.sum(column == max)) / float(numseq)  
    if percent_max >= 0.5:
        if max != 4:
            consensus_domain = max
    
    # Figure out where we are in the GPCR, based on consensus and previous columns
    if consensus_domain != GPCR_STRUCTURE[domain_index]:
        if consensus_domain == GPCR_STRUCTURE[domain_index + 1]:
            domain_index += 1
        elif consensus_domain != 3:
            raise AssertionError("\n\nWe have gotten lost in the GPCR structure. This is not ideal.")
    return domain_index, consensus_domain





def save_good_seqs(aln, keep_list, outfile):
    ''' Create resulting protein alignment for a given iteration based on keep_list'''
    count = 0
    outf = open(outfile, 'w')
    for record in aln:
        if keep_list[count]:
            newseq = str(record.seq).replace('-','')
            outf.write(">"+str(record.id)+"\n"+newseq+"\n")
        count += 1
    outf.close()

        

            
    
def save_nuc_struc_alns(prot_aln, nuc_infile, nuc_outfile, struc_aln, domain_outfile):
    ''' Use the final protein alignment file to create the final nucleotide alignment file. 
        For each nucleotide sequence, determine if we will keep it. If so, convert it to the aligned version, using the protein alignment.
        We also save, in fasta format, the protein alignment in domain form.
        Note that input prot_aln is an object.
    '''

    raw_nuc = list(SeqIO.parse(nuc_infile, "fasta"))
    outfile_nuc = open(nuc_outfile, "w")    
    outfile_dom = open(domain_outfile, "w")
    
    i = 0
    for prot in prot_aln:
        id = str(prot.id)
        prot_seq = str(prot.seq)
        struc_seq = "".join(map(str,struc_aln[i])).translate(TRANSLATOR)
        outfile_dom.write(">" + id + "\n" + struc_seq + "\n")
        for nuc in raw_nuc:
            if str(nuc.id) == id:
                old_nuc_seq = str(nuc.seq)
                new_nuc_seq = pal2nal_seq(prot_seq, old_nuc_seq)
                outfile_nuc.write(">" + id + "\n" + new_nuc_seq + "\n")
                break
        i += 1
    outfile_nuc.close()
    outfile_dom.close()


def pal2nal_seq(protseq, old_nucseq):
    ''' Convert a nucleotide sequence to aligned version, using the protseq (which is aligned!) as a guide for gap insertion. Works for a single pair of sequences, not a whole alignment.'''
    
    assert(len(old_nucseq) == len(protseq.replace("-",""))*3), "\n\nThese protein and nucleotide sequences have no business being compared."
    new_nucseq = ''
    start = 0
    end = 3
    for position in protseq:
        if position == '-' or position == '?':
            new_nucseq += 3*position
        elif position == 'B' or position == 'X' or position == 'Z':
            new_nucseq += 'NNN'
        else:
            new_nucseq += old_nucseq[start:end]
            start+=3
            end+=3
    assert(len(new_nucseq) == len(protseq)*3), "\n\nCould not convert nucleotide sequence to aligned version."
    return new_nucseq





def mask_alignment(numseq, numcol, seq_domains, prot_aln, nuc_alnfile, prot_mask_outfile, nuc_mask_outfile):
    ''' Mask protein residues which cannot be aligned structurally. This includes ambiguous!
        NOTE: prot_aln is alignment object and nuc_seqfile is unmasked nucleotide alignment file.
    '''
    
       
    struc_aln, num_aa = aln_to_gpcrhmm(prot_aln, numseq, numcol, seq_domains)
    domain_index = 0
    consensus_domain = 0
    full_structure = ''
    
    # Create a string representing the consensus domains
    for column in struc_aln.T:
        domain_index, consensus_domain = determine_consensus(numseq, column, domain_index, consensus_domain)
        full_structure += str(consensus_domain)
    assert(len(full_structure) == len(struc_aln[0])), "\n\nCould not compute the final consensus structure"

    # Now mask the residues based on agreement with consensus domains. If consensus domain is O,I,M, mask any violating residues. Elif consensus domain is A (ambiguous), mask the whole column. 
    nuc_aln = AlignIO.read(nuc_alnfile, "fasta")
    outf_p = open(prot_mask_outfile, "w")
    outf_n = open(nuc_mask_outfile, "w") 
    
    # Loop over sequences and whether to keep or mask each position, based on full_structure
    for i in range(numseq):
                  
        # Retrieve domains and protein and nucleotide unmasked sequences
        domains = struc_aln[i]
        prot_record = prot_aln[i]
        nuc_record = nuc_aln[i]
        prot_id = str(prot_record.id)
        nuc_id = str(nuc_record.id)
        assert(prot_id == nuc_id), "\n\nProtein and nucleotide records do not match in mask_alignment."
        prot_seq = str(prot_record.seq)
        nuc_seq = str(nuc_record.seq).replace("-", "") # remove gaps for pal2nal. not pretty, but I'm at peace with that.

        # For each entry in the column, determine whether to mask or not
        prot_mask = ''
        for j in range(numcol):
            if domains[j] != 4 and (domains[j] == 3 or domains[j] != int(full_structure[j])):
                prot_mask += '?'
            else:
                prot_mask += prot_seq[j]
        assert(len(prot_seq) == len(prot_mask)), "\n\nProtein sequence in the alignment not properly masked."
                
        # Convert to nucleotide mask
        nuc_mask = pal2nal_seq(prot_mask, nuc_seq)
        
        # Write masked to file
        outf_p.write(">" + prot_id + "\n" + prot_mask + "\n")
        outf_n.write(">" + nuc_id + "\n" + nuc_mask + "\n")             
    outf_p.close()
    outf_n.close()
    
    return full_structure



def save_naive_alignments(protaln, nucfile, outfile_protein, outfile_nuc):
    ''' Save first iteration mafft alignment. Naive! '''

    nuc_records = list(SeqIO.parse(nucfile, "fasta"))
    outf = open(outfile_nuc, "w")
    for i in range(len(nuc_records)):
        assert(str(protaln[i].id) == str(nuc_records[i].id)), "\n\nCan't make naive nucleotide alignment since protein and nucleotide ids don't match."
        newseq = pal2nal_seq( str(protaln[i].seq), str(nuc_records[i].seq) )
        outf.write(">" + str(protaln[i].id) + "\n" + newseq + "\n")
    outf.close()
    AlignIO.write(protaln, outfile_protein, "fasta")
    
    
    
    
    

def main():
    
    # I/O definitions
    assert(len(sys.argv) == 4), "\n\nUsage: python obtain_alignment.py <protein_input_file> <nucleotide_input_file> <gpcrhmm_directory>"
    infile_protein         = sys.argv[1]
    infile_nuc             = sys.argv[2]
    gpcrhmm_directory      = sys.argv[3]
    if gpcrhmm_directory[-1] != "/":
        gpcrhmm_directory += "/" 
    outfile_protein_naive  = "protein_aln_naive.fasta"
    outfile_nuc_naive      = "nucleotide_aln_naive.fasta"
    outfile_protein_masked = "protein_aln_masked.fasta"
    outfile_nuc_masked     = "nucleotide_aln_masked.fasta"
    outfile_protein        = "protein_aln.fasta"
    outfile_nuc            = "nucleotide_aln.fasta"
    outfile_domain_aln     = "domain_alignment.fasta"
    outfile_cons_domains   = "domain_consensus.txt" 
    mafft_out = "out.fasta"
    mafft_in  = "in.fasta"
    shutil.copy(infile_protein, mafft_in)
    
    # Before aligning, convert all sequences to their domains. As we iterate and remove sequences, we'll need to delete corresponding entries from this list to keep everything ordered in concert
    print "Converting sequences to domains"
    seq_domains = seq_to_gpcrhmm(infile_protein, gpcrhmm_directory)

    num_iterations = 1
    num_discarded = 1
    while num_discarded > 0:
        
        print
        print "Beginning iteration", num_iterations
        
        # Align with mafft
        print "Aligning"
        call_mafft(mafft_in, mafft_out, options='--quiet --thread 3')

        # Grab the output alignment to pass around to various functions
        prot_aln_obj = AlignIO.read(mafft_out, 'fasta')
        numseq = len(prot_aln_obj)
        numcol = len(prot_aln_obj[0])  

        # Check if we should save naive alignments (first iteration)
        if num_iterations == 1:
            save_naive_alignments(prot_aln_obj, infile_nuc, outfile_protein_naive, outfile_nuc_naive)

        # Convert protein alignment to gpcrhmm alignment. Returns protein alignment object, struc_aln numpy array of alignment
        print "Converting protein alignment to corresponding GPCR domains"
        struc_aln, num_aa = aln_to_gpcrhmm(prot_aln_obj, numseq, numcol, seq_domains)
       
        
        # Discard misaligned sequences and save unaligned protein sequences to file for next iteration
        print "Finding misaligned sequences to discard"
        keep_list = build_keep_list(struc_aln, numseq, num_aa)
        save_good_seqs(prot_aln_obj, keep_list, mafft_in)
        
        # Remove the structures records from seq_domains for those sequences which we are removing
        to_delete = list((keep_list == False).nonzero()[0])
        for offset, index in enumerate(to_delete):
            index -= offset
            del seq_domains[index]
            
        # Determine how many sequences were discarded
        num_discarded = np.sum(keep_list == False) # the number of sequences discarded is number of "False" entries in keep_list
        print num_discarded, "sequence(s) discarded"
        
        num_iterations += 1
    
    print "\nFinished protein alignment after", num_iterations-1, "iteration(s). Now saving final protein, nucleotide, and domain alignment files."
    shutil.copy(mafft_out, outfile_protein)
    save_nuc_struc_alns(prot_aln_obj, infile_nuc, outfile_nuc, struc_aln, outfile_domain_aln)
    
    print "\nCreating masked nucleotide and protein alignment files."
    full_structure = mask_alignment(numseq, numcol, seq_domains, prot_aln_obj, outfile_nuc, outfile_protein_masked, outfile_nuc_masked)
    outf = open(outfile_cons_domains, "w")
    outf.write(full_structure.translate(TRANSLATOR))
    outf.close()
    
    os.remove(mafft_in)
    os.remove(mafft_out)
    
        
main()
    
    
    
    
    
    
    
    
    
    
    

