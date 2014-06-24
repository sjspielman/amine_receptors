from Bio import SeqIO, AlignIO
import re
import sys
import os
import subprocess



gpcr_struc = 'OMIMOMIMOMIMOMIX' # need the X as as hack

########################## CONSENSUS STRUCTURE FUNCTIONS #################################
def maxColDict(column):
    ''' make dictionary to retrieve max entry in column.'''    
    dict = {}
    for entry in column:
        if entry in dict:
            dict[entry] += 1
        else:
            dict[entry] = 1
    mymax = max(dict.values())
    return dict, mymax
    
    
def findGPCRindex(column, index, c):
    ''' Where are we in the GPCR? '''
    dict, mymax = maxColDict(column)
    
    for entry in dict:
        if dict[entry] == mymax:
             maxstruc = entry
             break
    if maxstruc == gpcr_struc[index+1]:
        return index+1
    else:
        return index

def keepCol(col, domain, colindex):
    discard = []
    for r in range(len(col)):
        if col[r] == 'A' or col[r] == '-' or col[r] == domain:
            continue
        else:
            discard.append(r)
    return discard


def buildDiscardDict(seqs):
    ''' For each column, check to see if row shows proper domain.
        Build up dictionary of sequences which should be removed.
    '''
    index = 0
    full = ''
    discard_dict = {}
    for c in range(len(seqs[0])):
        col = ''
        for row in seqs:
            col += row[c]
        index = findGPCRindex(col, index, c)
        discard = keepCol(col, gpcr_struc[index], c)
        for entry in discard:
            if entry in discard_dict:
                discard_dict[entry] += 1
            else:
                discard_dict[entry] = 1
    return discard_dict
    
def saveGoodSeqs(seqfile, discard_dict, outfile, discardif):   
    numdisc = 0
    outf = open(outfile, 'w')
    raw = list(SeqIO.parse(seqfile, 'fasta'))
    for i in range(len(raw)):
        try:
            thresh = len(str(raw[i].seq).replace('-','')) * discardif
            if discard_dict[i] <= thresh:
                outf.write(">"+str(raw[i].id)+"\n"+str(raw[i].seq)+"\n")
            else:
                numdisc += 1
                
        except:
            outf.write(">"+str(raw[i].id)+"\n"+str(raw[i].seq)+"\n")
    outf.close()
    print "discarded",numdisc,"sequences"

    
def consStruc(strucfile, seqfile, outfile, discardif):
    
    # Grab seq2gpcrhmm "alignment"
    aln = list(SeqIO.parse(strucfile, 'fasta'))
    seqs = []
    for entry in aln:
        seqs.append(str(entry.seq))
    
    # Determine which sequences are structurally out of frame
    discard_dict = buildDiscardDict(seqs)
    
    # Save the good ones
    saveGoodSeqs(seqfile, discard_dict, outfile, discardif)

##########################################################################################



################################ SEQ2GPCRHMM FUNCTIONS ###################################


def probs2letter(inner, outer, membr, cutoff):
    ''' Turn probabilities going down the protein into a letter signifying domain.
        I = inner. O = outer. M = membrane. A = ambiguous.
        Return list of letters to replace alignment amino acids with.
    '''
    domain = []
    for i in range(len(inner)):
        temp = [ inner[i], outer[i], membr[i] ]
        if max(temp) >= cutoff:        
            which = temp.index(max(temp))
            if which == 0:
                domain.append('I')
            elif which == 1:
                domain.append('O')
            elif which == 2:
                domain.append('M')
            else:
                raise AssertionError("How in the hell did we get here?!!???!!!!")
        else:
            domain.append('A')
    assert(len(domain) == len(inner)), "A domain was not identified for all alignment positions. Womp womp."
    return domain            

def parseGPCRHMM(file, cutoff):
    ''' Parse a gpcrhmm file and return a list of domain letters to replace amino acid sequence with. ''' 
    hmm_file = open(file, 'rU')
    hmm_lines = hmm_file.readlines()
    hmm_file.close()
    
    innerList = []
    outerList = []
    membrList = []
    for line in hmm_lines:
        parseLine=re.search('^\d+\s\w\t(.+)\t(.+)\t(.+)\t(.+)$', line)
        if parseLine:
            innerList.append( float(parseLine.group(1)) )
            outerList.append( float(parseLine.group(2)) + float(parseLine.group(4)) )
            membrList.append( float(parseLine.group(3)) )
    assert(len(innerList) == len(outerList) == len(membrList) ), "GPCRHMM file not properly parsed."
    newLetters = probs2letter(innerList, outerList, membrList, cutoff)
    return newLetters    

def replaceSeq(aaseq, letters):
    ''' Replace the amino acid string (which includes gaps!!) with the domain letters. '''
    nogap_count = 0
    newseq = ''
    for position in aaseq:
        if position == '-':
            newseq += '-'
        else:
            newseq += letters[nogap_count]
            nogap_count += 1
    assert(len(newseq) == len(aaseq)), "AA sequence not properly replaced with domain indicators."
    return newseq

def amino2domain(infile, outfile, hmm_dir, cutoff):
    ''' Replace a GPCR alignment amino acids with domain.
        infile  = alignment infile
        outfile = outfile with new domain letters in lieu of amino acids
        hmm_dir = directory to the gpcrhmm files
        cutoff  = posterior probability threshold for calling domain.
    '''
    aln = AlignIO.read(infile, 'fasta')
    outaln = open(outfile, 'w')
    for record in aln:
        id = str(record.id)
        seq = str(record.seq)
        hmmFile = hmm_dir + id + '.txt'
        #print id
        newLetters = parseGPCRHMM(hmmFile, cutoff)
        finalSeq = replaceSeq(seq, newLetters)
        outaln.write(">"+id+"\n"+finalSeq+"\n")
    outaln.close()




