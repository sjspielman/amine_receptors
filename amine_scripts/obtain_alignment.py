# SJS 9/6/14 - 9/7/14.

import sys
import re
import os
import subprocess
import shutil
import numpy as np
from Bio import AlignIO, SeqIO

np.set_printoptions(threshold='nan')

################### GLOBALS ####################
GPCR_DOMAIN_DICT   = {'O':0, 'I':1, 'M':2, 'A':3} # note that A means ambiguous. Also note that -1 codes for a gap.
GPCR_STRUCTURE     = '0212021202120213' # we use the 3 as a "hack" to terminate the protein.
PERCENT_MISALIGNED = 0.05  # We discard sequences where >=5% of positions do not match consensus structure
GPCRHMM_THRESHOLD  = 0.5   # We assign domains if posterior probability is >= 0.5. Otherwise, we call it ambiguous.
GPCRHMM_DIRECTORY  = "/Users/sjspielman/Dropbox/GPCRs/Amine/gpcrhmm_records/"
################################################


def call_mafft(infile, outfile, options=''):
    ''' Infile contains unaligned protein fasta sequences to be sent to outfile. Options can be given to mafft as well.'''
    mafft = "mafft " + options + " " + infile + " > " + outfile
    runmafft = subprocess.call(mafft, shell=True)
    assert(runmafft == 0), "Mafft failed to run."

    
def aln_to_gpcrhmm(alnfile):
    ''' Convert a protein alignment file to a "structural" alignment. Replace amino acid residues with domain (O, I, M, A).'''
    aln = AlignIO.read(alnfile, 'fasta')    
    struc_aln = np.zeros( [len(aln), len(aln[0])], dtype = 'int8')
    i=0
    for record in aln:
        id = str(record.id)
        seq = str(record.seq)
        num_aa = len(seq) - seq.count('-')
        gpcrhmm_file = GPCRHMM_DIRECTORY +  id.split("_")[0]+"_" + id.split("_")[1] + '.txt'
        domains = parse_gpcrhmm(gpcrhmm_file, num_aa)
        domain_sequence = replace_residues(seq, domains)
        struc_aln[i] = domain_sequence
        i += 1
    return aln, struc_aln


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
            domains += '-1'
    assert(len(domains) == length), "A domain was not identified for all alignment positions. Womp womp."
    return domains
   
   
def replace_residues(seq, domains):
    ''' Replace the amino acid string (which includes gaps!!) with the domain letters. '''
    nogap_count = 0
    total_count = 0
    newseq = np.zeros(len(seq), dtype = 'int8')
    newseq[newseq == 0] = 5 
    
    for position in seq:
        if position == '-':
            newseq[total_count] = -1
        else:
            newseq[total_count] = GPCR_DOMAIN_DICT[ domains[nogap_count] ]
            nogap_count += 1
        total_count += 1
    assert(newseq.all() != 5), "Could not convert sequence to domains."
    return newseq




def build_keep_list(seq_array):
    
    domain_index = 0
    keep_list = np.zeros(len(seq_array)) # Each element is a row in the alignment. The resulting numbers represent how many columns in that sequence violate consensus.
    for column in seq_array.T:
        
        # Determine the column consensus. Note that we should not consider gaps when determining column consensus!
        unique, positions = np.unique(column[column != -1], return_inverse=True)
        consensus_domain = unique[np.bincount(positions).argmax()]
        
        # Figure out where we are in the GPCR, based on consensus and previous columns
        if consensus_domain == GPCR_STRUCTURE[domain_index+1]:
            domain_index += 1

        # Add 1 to keep_list for those rows which do not match consensus
        masks = [column != -1, column != 3, column != consensus_domain] 
        keep_list += reduce(np.logical_and, masks)
    
    # Return boolean array for whether we should keep each sequence.
    keep_list /= float(len(seq_array[0]))
    return keep_list < PERCENT_MISALIGNED

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
            
    
def save_nucleotide_aln(protein_file, nuc_infile, nuc_outfile):
    ''' Use the final protein alignment file to create the final nucleotide alignment file. 
        For each nucleotide sequence, determine if we will keep it. If so, convert it to the aligned version, using the protein alignment.
    '''

    prot_aln = AlignIO.read(protein_file, "fasta")
    raw_nuc = list(SeqIO.parse(nuc_infile, "fasta"))
    outfile = open(nuc_outfile, "w")    

    for prot in prot_aln:
        id = str(prot.id)
        prot_seq = str(prot.seq)
        for nuc in raw_nuc:
            if str(nuc.id) == id:
                old_nuc_seq = str(nuc.seq)
                new_nuc_seq = ''
                start = 0
                end = 3
                for position in prot_seq:
                    if position == '-':
                        new_nuc_seq += '---'
                    elif position == 'B' or position == 'X' or position == 'Z':
                        new_nuc_seq += 'NNN'
                    else:
                        new_nuc_seq += old_nuc_seq[start:end]
                        start+=3
                        end+=3
                assert(len(new_nuc_seq) == len(prot_seq)*3), "Could not convert nucleotide sequence to aligned version."
                outfile.write(">" + id + "\n" + new_nuc_seq + "\n")
                break
    outfile.close()
                
 
    
def main():
    
    assert(len(sys.argv) == 5), "\n\nUsage: python obtain_alignment.py <infile_protein> <infile_nucleotide> <outfile_protein> <outfile_nucleotide>."
    infile_protein = sys.argv[1]
    infile_nucleotide = sys.argv[2]
    outfile_protein = sys.argv[3]
    outfile_nucleotide = sys.argv[4]
    
    mafft_out = "out.fasta"
    mafft_in  = "in.fasta"
    shutil.copy(infile_protein, mafft_in)
    
    count = 0
    num_discarded = 1
    while num_discarded > 0:
        
        print
        print "Beginning iteration", count
        
        # Align with mafft
        print "Aligning"
        call_mafft(mafft_in, mafft_out, options='--quiet')

        # Convert protein alignment to gpcrhmm alignment. Returns protein alignment object, struc_aln numpy array of alignment
        print "Converting protein alignment to corresponding GPCR domains"
        prot_aln, struc_aln = aln_to_gpcrhmm(mafft_out)

        # Discard misaligned sequences and save unaligned protein sequences to file for next iteration
        print "Finding misaligned sequences to disard"
        keep_list = build_keep_list(struc_aln)
        save_good_seqs(prot_aln, keep_list, mafft_in)
    
        # Determine how many sequences were discarded
        num_discarded = np.sum(keep_list == False) # the number of sequences discarded is number of "False" entries in keep_list
        print num_discarded, "sequence(s) discarded"
        
        count += 1
    
    print "\nFinished protein alignment after", count-1, "iteration(s). Now saving final alignment files."
    # After final sequence set has been determined, create the final protein and nucleotide alignment files
    shutil.copy(mafft_out, outfile_protein)
    save_nucleotide_aln(outfile_protein, infile_nucleotide, outfile_nucleotide)
    
    os.remove(mafft_in)
    os.remove(mafft_out)
    
        
main()
    
    
    
    
    
    
    
    
    
    
    

