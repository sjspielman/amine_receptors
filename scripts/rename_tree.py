# SJS 10/1/14-10/2/14. Renames taxa in a given tree for easier post-processesing (as in, manual coloring) for figure creation.
# The values in dictionary named "final_amine_names" can be used to search/find taxa and color accordingly in FigTree or similar program.
# USAGE python rename_tree.py <input_tree> <output_tree> 

import os
import re
import sys
from Bio import SeqIO
base_directory = "/Users/sjspielman/Dropbox/GPCRs/Amine/"
ncbi_directory = base_directory + "ncbi_records/"
tree_directory = base_directory + "phylogenies/"


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
                 
                


def rename_taxon(tree_string, index):
    find = re.search('^([XN]P_\d+\.\d)(_[XN]M_\d+\.\d)', tree_string)
    assert(find is not None), "Can't parse taxon name."
    protid = find.group(1)
    fullid = find.group(1) + find.group(2)
    
    ncbi = SeqIO.read(ncbi_directory + protid + ".txt", "gb")
    desc = "fudgetext " + str(ncbi.description).upper().replace('(', '').replace(')','')
    subtype = 'unknown'; key = 'unknown';
    for type in amine_subtypes:
        if type in desc:
            try:
                key = type
                subtype = amine_subtypes[type].match(desc).group(1) # Validated that works for all sequences which are not unknowns
                break
            except:
                pass
    new_name = fullid + '_' + final_amine_names[key] + '_' + subtype
    new_index = index + len(fullid)
    return new_name, new_index
    

def rename_tree(infile, outfile):
    inf = open(infile, 'r')
    tree_string = inf.read()
    inf.close()
    tree_string = re.sub(r"\s", "", tree_string)
    treelen = len(tree_string)
    out_tree = ''
    
    index = 0
    while index < treelen:
        print index  #, out_tree
        if tree_string[index] == 'X' or tree_string[index] == 'N':
            add_string, index = rename_taxon(tree_string[index:], index)
            out_tree += add_string
        else:
            out_tree += tree_string[index]
            index += 1    
    
    outf = open(outfile, 'w')
    outf.write(out_tree)
    outf.close()
            
            
def main():
    assert(len(sys.argv) == 3), "Usage: python rename_tree.py <input_tree> <output_tree> . \nTree files are *assumed* to be in directory ~/Dropbox/Amine/GPCRs/phylogenies/ . "
    in_tree_file  = tree_directory + sys.argv[1]
    out_tree_file = tree_directory + sys.argv[2]
    assert(os.path.exists(in_tree_file)), "Your input tree file does not exist in ~/Dropbox/Amine/GPCRs/phylogenies/. Quitting."      
    rename_tree(in_tree_file, out_tree_file)
            
            
            
            