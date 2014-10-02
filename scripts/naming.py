# SJS 10/1/14. So far, determines the type, subtype for all sequences. Note that 13 of them are *incorrectly* annotated (come out as unknown_unknown).
# Returns, for each sequence, something like "5ht_1A". We want to use this script once trees are finished in order to color in figtree.

import re
from Bio import SeqIO
amine_subtypes = {'DOPAMINE':             re.compile('.+(D\d|\dD).+'),
                  'HISTAMINE':            re.compile('.+H(\d).+'), 
                  'TRACE':                re.compile('.+RECEPTOR ([1-9]\w*).+'),
                  'ADREN':                re.compile('.+(\w+[- ][1-9]\w*).+'),
                  'ACETYLC':              re.compile('.+RECEPTOR (M\d\w*).+'),
                  'MUSCARINIC':           re.compile('.+\s*MUSCARINIC[ A-Z]*(\d\w*).+'),
                  'SEROTONIN 5-HT':       re.compile('.+SEROTONIN 5-HT(\d).+'),
                  'SEROTONIN':            re.compile('.+RECEPTOR ([1-9]\w*).+'),
                  '5-HT':                 re.compile('.+5-HT(\d).+'), 
                  '5-HYDROXYTRYPTAMINE':  re.compile('.+RECEPTOR ([1-9]\w*).+'),
                  }
final_amine_names = {'DOPAMINE':            'dopa',
                     'HISTAMINE':           'hrh',  
                     'TRACE':               'trace',   
                     'ADREN':               'adr', 
                     'ACETYLC':             'cholin', 
                     'MUSCARINIC':          'cholin', 
                     'SEROTONIN 5-HT':      '5ht',
                     'SEROTONIN':           '5ht',
                     '5-HT':                '5ht',
                     '5-HYDROXYTRYPTAMINE': '5ht',
                     'unknown':             'unknown'
                    }
                 
                
                 
base_directory = "/Users/sjspielman/Dropbox/GPCRs/Amine/"
ncbi_directory = base_directory + "ncbi_records/"
aln_directory  = base_directory + "alignments/"
aln = list(SeqIO.parse(aln_directory + "protein_aln.fasta", "fasta"))

for record in aln:
    id = str(record.id)
    protid = id.split('_')[0] + '_' + id.split('_')[1] 
    
    ncbi = SeqIO.read(ncbi_directory + protid + ".txt", "gb")
    desc = "fudgetext " + str(ncbi.description).upper().replace('(', '').replace(')','')
    
    # Determine receptor type
    subtype = 'unknown'
    key = 'unknown'
    for type in amine_subtypes:
        if type in desc:
            try:
                key = type
                subtype = amine_subtypes[type].match(desc).group(1) # Validated that works for all sequences which are not unknowns
                break
            except:
                pass
    final_name = final_amine_names[key] + '_' + subtype
   
            
            
            
            
            
            
            
            
            
            
            