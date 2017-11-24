import numpy as np
from filenames import tree_newick, tree_sequence_alignment
from util import (read_tree_meta_data, write_tree_meta_data,
                  get_genes_and_alignments, generic_argparse)
from Bio import Phylo

def get_amino_acid_mutations(tree, fname):
    seqs = {}
    from Bio import SeqIO
    for seq in SeqIO.parse(fname, 'fasta'):
        seqs[seq.name] = seq

    muts = {}
    muts[T.root.name]=''
    for node in T.get_nonterminals():
        pseq = seqs[node.name]
        for c in node:
            cseq = seqs[c.name]
            muts[c.name]=','.join([anc+str(pos+1)+der
                        for pos, (anc, der) in enumerate(zip(pseq, cseq))
                        if anc!=der])

    return muts


if __name__ == '__main__':
    parser = generic_argparse("Assign amino acid mutations to the tree")
    args = parser.parse_args()
    path = args.path

    tree_meta = read_tree_meta_data(path)
    T = Phylo.read(tree_newick(path), 'newick')

    for gene, aln_fname in get_genes_and_alignments(path, tree=True):
        if gene!='nuc':
            muts = get_amino_acid_mutations(T, aln_fname)

        for node_name in tree_meta:
            tree_meta[node_name][gene+'_mutations'] = muts[node_name]
    write_tree_meta_data(path, tree_meta)
