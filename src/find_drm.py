#EBH 15 Feb 2018
#Identifies DRMs in VCF-files and writes them into tree_meta file.
#Also creates colours for them and write them to a new drm_color file which is
#later included in export_to_auspice if it exists

import numpy as np
import colorsys
from util import read_tree_meta_data, write_tree_meta_data, read_in_vcf, generic_argparse
from filenames import tree_vcf_alignment, ref_fasta, tree_meta_file_name, drm_color_maps


def read_in_DRMs(drm_file):
    import pandas as pd

    DRMs = {}
    drmPositions = []

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = m.GENOMIC_POSITION-1 #put in python numbering
        drmPositions.append(pos)

        if pos in DRMs:
            DRMs[pos]['base'].append(m.ALT_BASE)
            DRMs[pos]['AA'].append(m.SUBSTITUTION)
        else:
            DRMs[pos] = {}
            DRMs[pos]['base'] = [m.ALT_BASE]
            DRMs[pos]['drug'] = m.DRUG
            DRMs[pos]['AA'] = [m.SUBSTITUTION]
            DRMs[pos]['gene'] = m.GENE

    drmPositions = np.array(drmPositions)
    drmPositions = np.unique(drmPositions)
    drmPositions = np.sort(drmPositions)

    DRM_info = {'DRMs': DRMs,
            'drmPositions': drmPositions}

    return DRM_info


def drugTranslate(x):
    return {
        'RIF': 'Rifampicin',
        'FQ': 'Fluoroquinolones',
        'ETH': 'Ethionamide',
        'EMB': 'Ethambutol',
        'SM': 'Streptomycin',
        'PZA': 'Pyrazinamide',
        'CAP': 'Capreomycin',
        'INH': 'Isoniazid',
        'KAN': 'Kanamycin',
        'AK': 'Amikacin'
    }[x]


def get_N_HexCol(n=5):
    #modified from kquinn's answer on stackoverflow:
    #https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
    import colorsys

    HSV_tuples = [(x*1.0/n, 0.7, 0.7) for x in xrange(n)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append("#"+"".join(map(lambda x: chr(x).encode('hex'),rgb)))
    return hex_out


def find_drms(DRM_info, sequences):
    drmPositions = DRM_info['drmPositions']
    DRMs = DRM_info['DRMs']

    seqDRM = {}
    for key in drmPositions:
        for seq,v in sequences.iteritems():
            try:
                base = sequences[seq][key]
                if base in DRMs[key]['base']:
                    if seq not in seqDRM:
                        seqDRM[seq] = {}

                    i=0
                    while base != DRMs[key]['base'][i]:
                        i+=1
                    seqDRM[seq][DRMs[key]['AA'][i]] = DRMs[key]['drug']

            except KeyError, e:
                continue

    return seqDRM


def remove_old_DRM(tree_met):
    tree_met.pop('Rifampicin', None)
    tree_met.pop('Fluoroquinolones', None)
    tree_met.pop('Ethionamide', None)
    tree_met.pop('Ethambutol', None)
    tree_met.pop('Streptomycin', None)
    tree_met.pop('Pyrazinamide', None)
    tree_met.pop('Capreomycin', None)
    tree_met.pop('Isoniazid', None)
    tree_met.pop('Kanamycin', None)
    tree_met.pop('Amikacin', None)
    tree_met.pop('Drug_Resistance', None)


def add_drm_tree_meta(path, seqDRM):
    tree_meta = read_tree_meta_data(path)

    #add drug resistance to tree_meta, & make list for colouring
    drugMuts = {}
    drugMuts["Drug_Resistance"] = ['0']
    for seq,v in seqDRM.iteritems():
        #in case re-running, don't add mutations to old ones!
        remove_old_DRM(tree_meta[seq])
        tempList = {}
        for mut,drug in v.iteritems():
            drugs = drug.split(';')
            for drug in drugs:
                trDrug = drugTranslate(drug)
                if trDrug in tree_meta[seq]:
                    tree_meta[seq][trDrug] = ",".join([tree_meta[seq][trDrug],mut])
                else:
                    tree_meta[seq][trDrug] = mut

                if trDrug in drugMuts:
                    if tree_meta[seq][trDrug] not in drugMuts[trDrug]:
                        drugMuts[trDrug].append(tree_meta[seq][trDrug])
                else:
                    drugMuts[trDrug] = [ tree_meta[seq][trDrug] ]

                tempList[trDrug] = ""

        numResist = str(len(tempList))
        tree_meta[seq]["Drug_Resistance"] = numResist
        if numResist not in drugMuts["Drug_Resistance"]:
            drugMuts["Drug_Resistance"].append(numResist)

    #for any with no resistance, add a 0 to tree_meta
    for seq,v in tree_meta.iteritems():
        if 'Drug_Resistance' not in tree_meta[seq]:
            tree_meta[seq]["Drug_Resistance"] = '0'

    write_tree_meta_data(path, tree_meta)

    return drugMuts


if __name__ == '__main__':
    parser =  generic_argparse("Find drug resistance mutations according to supplied file. ONLY WORKS FOR VCF FILES.")
    parser.add_argument('--drm', type=str,
                        help="file of DRMs to find")

    args = parser.parse_args()
    path = args.path

    compress_seq = read_in_vcf(tree_vcf_alignment(path), ref_fasta(path), compressed=False)

    sequences = compress_seq['sequences']
    positions = compress_seq['positions']
    ref = compress_seq['reference']

    DRM_info = read_in_DRMs(args.drm)

    #Find the DRMs and store them
    seqDRM = find_drms(DRM_info, sequences)

    #Add DRMs to the tree_meta file, and also collect all unique options/combos sort
    #we can colour them
    drugMuts = add_drm_tree_meta(path, seqDRM)

    #Get colours and write them to a new color file that will be checked
    #in export_to_auspice and included in meta.json if present
    newColAdds = []
    drugMuts["Drug_Resistance"].sort()
    for drug,muts in drugMuts.iteritems():
        cols = get_N_HexCol(len(muts))
        drugCol = [ "\t".join([drug, muts[i], cols[i]]) for i in xrange(len(muts)) ]
        drugCol.append("")
        newColAdds = newColAdds + drugCol

    with open(drm_color_maps(path), 'w') as the_file:
        the_file.write("\n".join(newColAdds))


