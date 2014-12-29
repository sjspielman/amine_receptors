# SJS 10/1/14-10/2/14. Renames taxa in a given tree for easier post-processesing (as in, manual coloring) for figure creation.
# The values in dictionary named "final_amine_names" can be used to search/find taxa and color accordingly in FigTree or similar program.
# USAGE python rename_tree.py <input_tree> <output_tree> 

import os
import re
import sys
from Bio import SeqIO
base_directory = "/Users/sjspielman/Research/amine_receptors/analysis/"
ncbi_directory = base_directory + "ncbi_records/"
tree_directory = base_directory + "phylogeny/"


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
final_amine_names = {'DOPAMINE':            'drd',
                     'HISTAMINE':           'hrh',  
                     'TRACE':               'taar',   
                     'ADREN':               'adr', 
                     'ACETYLC':             'machr', 
                     'MUSCARINIC':          'machr', 
                     'SEROTONIN 5-HT':      '5ht',
                     'SEROTONIN':           '5ht',
                     '5-HT':                '5ht',
                     '5-HYDROXYTRYPTAMINE': '5ht',
                     'unknown':             'unknown'
                    }
                 
                 
misannotated = {"XP_005797918.1_XM_005797861.1": "drd-3", "XP_003967971.1_XM_003967922.1": "drd-3", "NP_001266433.1_NM_001279504.1": "machr-2", "XP_001520508.2_XM_001520458.3": "hrh-4", "XP_005282846.1_XM_005282789.1": "hrh-3", "XP_001920844.1_XM_001920809.1": "taar-12", "NP_001076571.1_NM_001083102.1": "taar-13", "XP_006014096.1_XM_006014034.1": "taar-4", "XP_003201718.2_XM_003201670.2": "taar-10", "NP_001076546.1_NM_001083077.1": "taar-10", "NP_001083418.1_NM_001089949.1": "adr-b", "NP_001103208.1_NM_001109738.1": "hrh-2", "NP_001124143.1_NM_001130671.1": "taar-12", "XP_001337671.1_XM_001337635.2": "taar-star", "XP_003976403.1_XM_003976354.1": "taar-star", "XP_005810466.1_XM_005810409.1": "taar-star", "XP_003454279.1_XM_003454231.1": "taar-star", "XP_004549625.1_XM_004549568.1": "taar-star", "XP_002935532.2_XM_002935486.2": "taar-star", "XP_006013317.1_XM_006013255.1": "taar-star", "XP_005510029.1_XM_005509972.1": "unknown", "XP_002187301.2_XM_002187265.2": "unknown", "XP_002937327.2_XM_002937281.2": "unknown", "XP_005045681.1_XM_005045624.1": "unknown", "XP_005144673.1_XM_005144616.1": "unknown", "XP_005229932.1_XM_005229875.1": "unknown", "XP_005428400.1_XM_005428343.1": "unknown", "XP_005490920.1_XM_005490863.1": "unknown", "XP_005518128.1_XM_005518071.1": "unknown", "XP_006111669.1_XM_006111607.1": "unknown", "XP_420867.2_XM_420867.4": "unknown"}               

def rename_taxon(tree_string, index, update):
    find = re.search('^([XN]P_\d+\.\d)(_[XN]M_\d+\.\d)', tree_string)
    assert(find is not None), "Can't parse taxon name."
    protid = find.group(1)
    fullid = find.group(1) + find.group(2)

    ncbi = SeqIO.read(ncbi_directory + protid + ".txt", "gb")
    tax = "".join(ncbi.annotations["taxonomy"]).replace("EukaryotaMetazoaChordataCraniataVertebrataEuteleostomi","")
    desc = "fudgetext " + str(ncbi.description).upper().replace('(', '').replace(')','')
    
    if update and fullid in misannotated:
        new_name = fullid + '_' + misannotated[fullid] + '_' + tax
    else:
        subtype = 'unknown'; key = 'unknown';
        for type in amine_subtypes:
            if type in desc:
                try:
                    key = type
                    subtype = amine_subtypes[type].match(desc).group(1) # Validated that works for all sequences which are not unknowns
                    break
                except:
                    pass
        subtype = subtype.replace('_','-').replace(' ','')
        new_name = fullid + '_' + final_amine_names[key] + '-' + subtype + '_' + tax

    new_index = index + len(fullid)
    return new_name, new_index
    

def rename_tree(infile, outfile, update):
    inf = open(infile, 'r')
    tree_string = inf.read()
    inf.close()
    tree_string = re.sub(r"\s", "", tree_string)
    treelen = len(tree_string)
    out_tree = ''
    
    index = 0
    while index < treelen:
        print index
        #print index  #, out_tree
        if tree_string[index] == 'X' or tree_string[index] == 'N':
            add_string, index = rename_taxon(tree_string[index:], index, update)
            out_tree += add_string
        else:
            out_tree += tree_string[index]
            index += 1    
    outf = open(outfile, 'w')
    outf.write(out_tree)
    outf.close()
    
            
            
def main():
    assert(len(sys.argv) == 3 or len(sys.argv) == 4), "Usage: python rename_tree.py <input_tree> <output_tree> <update_misannotated>. \nTree files are *assumed* to be in directory ../analysis/phylogeny/. Final argument (True/False) is optional (default False). If true, output tree will have updated annotations as described in paper. "
    in_tree_file  = tree_directory + sys.argv[1]
    out_tree_file = tree_directory + sys.argv[2]
    try:
        update = bool(sys.argv[3])
    except:
        update = False
    assert(os.path.exists(in_tree_file)), "Your input tree file does not exist. Quitting."      
    rename_tree(in_tree_file, out_tree_file, update)

main()
            
            
            