import numpy as np
import os, shutil, glob
from util import generic_argparse, load_features, read_tree_meta_data, write_tree_meta_data, get_genes_and_alignments, write_json
from translate import translate_feature, get_amino_acid_mutations, assign_clades
from export_to_auspice import attach_tree_meta_data, tree_to_json, tree_layout
from filenames import ref_alignment, tree_newick, tree_sequence_alignment, tree_json
from Bio import Phylo, AlignIO


if __name__ == '__main__':
    parser = generic_argparse("Infer mutations for one segment on the tree of another")
    parser.add_argument('--prefix', required=True,
                        help="prefix for json files that are passed on to auspice (e.g., zika.fasta)")

    args = parser.parse_args()
    path = args.path
    base_path = '/'.join(path.rstrip('/').split('/')[:-1]) + '/'

    segments = map(lambda x:x.split('/')[-1], filter(lambda x:os.path.isdir(x) and 'meta' not in x, glob.glob(base_path+'/*')))

    # load HA tree meta data with clade
    ha_tree_meta = assign_clades(path)
    for segment in segments:
        tree_meta = read_tree_meta_data(base_path+segment)
        for node in tree_meta:
            if node.startswith("NODE"):
                continue
            elif node in ha_tree_meta and ("named_clades" in ha_tree_meta[node]):
                tree_meta[node]["named_clades"] = ha_tree_meta[node]["named_clades"]

        write_tree_meta_data(base_path+segment, tree_meta)