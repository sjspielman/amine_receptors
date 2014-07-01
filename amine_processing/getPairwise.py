## SJS 5/13/14. Script to calculate pairwise similarities for a sequence file.
## Build a temp full alignment with linsi (assumes we have <=200 seqs, which we probably do)
## Build a list of sequences that need to be removed based on a sequence similarity cutoff, and ultimately remove those sequences.


import sys
import re
import subprocess
from Bio import SeqIO


def Align(infile, outfile):
    alignCall = 'linsi ' + infile + ' > ' + outfile
    check = subprocess.call(alignCall, shell=True)
    assert(check == 0), "Alignment call did not work."
    
    
    
def calcPairSim(seq1, seq2, id1, id2, cutoff, cullList):
    ''' Calculate pairwise similarity between two sequences. Cull one if necessary.'''
    same = 0.
    shared_length = 0.
    for i in range(len(seq1)):
        if seq1[i] == '-' and seq2[i] == '-':
            continue
        else:
            shared_length += 1.
            if seq1[i] == seq2[i]:
                same += 1.
    sim = same/shared_length  
    
    if sim >= cutoff:
        if id1 not in cullList and id2 not in cullList:
            if "NP" in id1 and "NP" not in id2:
                cullList.append(id2)
            elif "NP" in id2 and "NP" not in id1:
                cullList.append(id1)
            elif "unknown" in id1 and "unknown" not in id2:
                cullList.append(id1)
            elif "unknown" in id2 and "unknown" not in id1:
                cullList.append(id2)
            else:
                len_seq1 = len(seq1) - seq1.count('-')
                len_seq2 = len(seq2) - seq2.count('-')
                if len_seq1 >= len_seq2:
                    cullList.append(id2)
                else:
                    cullList.append(id1)
    return cullList
    
    
def main(seqfile, cutoff):

    alnfile = 'amine_aa.aln'
    #Align(seqfile, alnfile)
    allSeqs = list(SeqIO.parse(alnfile, 'fasta'))
    
    cullList = []
    
    # create cullList to figure out which seqs to remove based on similarity
    for n in range(len(allSeqs)):
        print n
        print cullList
        entry1 = allSeqs[n]
        id1 = str(entry1.id)
        seq1 = str(entry1.seq)

        for m in range(n+1, len(allSeqs)):
            if id1 in cullList:
                break
            entry2 = allSeqs[m]
            if entry1.id != entry2.id and entry2.id not in cullList:
                id2 = str(entry2.id)
                seq2 = str(entry2.seq)
                cullList = calcPairSim(seq1, seq2, id1, id2, cutoff, cullList)
                    
    # save sequences we want to keep
    outf = open(outfile, 'w')
    for entry in allSeqs:
        if str(entry.id) not in cullList:
            seq = str(entry.seq).translate(None, '-') # rm gaps
            outf.write(">"+str(entry.id)+"\n"+str(seq)+"\n")
    outf.close()
    return cullList
    
    
assert(len(sys.argv) == 4), "Usage: python getPairwise.py <infile> <outfile> <cutoff>. Outfile is a fasta file containing sequences we are keeping. Cutoff is decimal cutoff for similarity. Above cutoff will not be saved."
infile = sys.argv[1]
outfile = sys.argv[2]
cutoff = float(sys.argv[3])
cullList = main(infile, cutoff)

print "Number of sequences removed due to similarity:", len(cullList)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    