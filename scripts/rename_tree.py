# SJS 10/1/14-10/2/14. Renames taxa in a given tree for easier post-processesing (as in, manual coloring) for figure creation.
# The values in dictionary named "final_amine_names" can be used to search/find taxa and color accordingly in FigTree or similar program.
# USAGE python rename_tree.py <input_tree> <output_tree> <update_taxonomy>    (last arg is optional bool)

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
final_amine_names = {'DOPAMINE':            'DRD',
                     'HISTAMINE':           'HRH',  
                     'TRACE':               'TAAR',   
                     'ADREN':               'ADR', 
                     'ACETYLC':             'CHRM', 
                     'MUSCARINIC':          'CHRM', 
                     'SEROTONIN 5-HT':      'HTR',
                     'SEROTONIN':           'HTR',
                     '5-HT':                'HTR',
                     '5-HYDROXYTRYPTAMINE': 'HTR'} 

                 
                 
misannotated = {"XP_005797918.1_XM_005797861.1": "DRD3", "XP_003967971.1_XM_003967922.1": "DRD3", "NP_001266433.1_NM_001279504.1": "MACHR2", "XP_001520508.2_XM_001520458.3": "HRH4", "XP_005282846.1_XM_005282789.1": "HRH3", "XP_001920844.1_XM_001920809.1": "TAAR12", "NP_001076571.1_NM_001083102.1": "TAAR13", "XP_006014096.1_XM_006014034.1": "TAAR4", "XP_003201718.2_XM_003201670.2": "TAAR10", "NP_001076546.1_NM_001083077.1": "TAAR10", "NP_001083418.1_NM_001089949.1": "ADRB", "NP_001103208.1_NM_001109738.1": "HRH2", "NP_001124143.1_NM_001130671.1": "TAAR12", "XP_001337671.1_XM_001337635.2": "TAARstar", "XP_003976403.1_XM_003976354.1": "TAARstar", "XP_005810466.1_XM_005810409.1": "TAARstar", "XP_003454279.1_XM_003454231.1": "TAARstar", "XP_004549625.1_XM_004549568.1": "TAARstar", "XP_002935532.2_XM_002935486.2": "TAARstar", "XP_006013317.1_XM_006013255.1": "TAARstar", "XP_005510029.1_XM_005509972.1": "unknown", "XP_002187301.2_XM_002187265.2": "unknown", "XP_002937327.2_XM_002937281.2": "unknown", "XP_005045681.1_XM_005045624.1": "unknown", "XP_005144673.1_XM_005144616.1": "unknown", "XP_005229932.1_XM_005229875.1": "unknown", "XP_005428400.1_XM_005428343.1": "unknown", "XP_005490920.1_XM_005490863.1": "unknown", "XP_005518128.1_XM_005518071.1": "unknown", "XP_006111669.1_XM_006111607.1": "unknown", "XP_420867.2_XM_420867.4": "unknown"}               

def rename_taxon(tree_string, index, update):
    find = re.search('^([XN]P_\d+\.\d)(_[XN]M_\d+\.\d)', tree_string)
    assert(find is not None), "Can't parse taxon name."
    protid = find.group(1)
    fullid = find.group(1) + find.group(2)
    type = ''
    ncbi = SeqIO.read(ncbi_directory + protid + ".txt", "gb")
    
    if update and fullid in misannotated:
        new_name = fullid + '-' + misannotated[fullid]

    else:
        for feat in ncbi.features:
            if feat.type == 'CDS':
                type = feat.qualifiers['gene'][0].upper()

        # grab info from description if gene poorly defined       
        if "LOC" in type or "SI" in type:
            desc = "fudgetext " + str(ncbi.description).upper().replace('(', '').replace(')','')
            subtype = 'unknown'; key = 'unknown';
            for amine in amine_subtypes: 
                if amine in desc:
                    try:
                        key = amine
                        subtype = amine_subtypes[amine].match(desc).group(1) # Validated that works for all sequences which are not unknowns
                        break
                    except:
                        pass
            type = final_amine_names[key]+subtype
        
        # adrenergic receptors are noncanonically annotated in their description, so grab again from protein feature
        if "ADR" in type:
            for feat in ncbi.features:
                if feat.type == "Protein":
                    temp = feat.qualifiers['product'][0].upper()
                    if "BETA" in temp:
                         type = 'ADRB'
                    else:
                        find_adra = re.search('ALPHA[- ](\d)[A-Z]', temp)
                        if find_adra:
                            type = 'ADRA'+find_adra.group(1)

    
        
        ####### type cleanups ########
        # remove symbols
        type = type.replace('_','').replace('-','')

        # drd-specific
        type = type.replace('DRDD', 'DRD')
        
         # extra trailing letters
        find = re.search("^(.+\d)[A-Z]+$", type)
        if find:
            type = find.group(1)
        
        # some htr names have these letters which I want gone.
        find = re.search("(.+)CL\d+$", type)
        if find:
            type = find.group(1)
        
        # one random adr looks this way
        if type == 'ADB4C':
            type = 'ADRB'
        ########################################
        
        new_name = fullid + '_' + type

    # THIS LINE WILL BASICALLY WRITE THE SEQUENCE_DESCRIPTIONS.TXT FILE !!!
    print fullid + '\t' + type + '\t' + ".".join(ncbi.annotations["taxonomy"])#.replace("EukaryotaMetazoaChordataCraniataVertebrataEuteleostomi","")
    
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
        #print index
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
            
            
            