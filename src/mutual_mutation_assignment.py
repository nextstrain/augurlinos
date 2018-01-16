import numpy as np
import os, shutil, glob
from util import generic_argparse, load_features, read_tree_meta_data, write_tree_meta_data, get_genes_and_alignments, write_json
from translate import translate_feature, get_amino_acid_mutations, assign_clades
from export_to_auspice import attach_tree_meta_data, tree_to_json, tree_layout
from filenames import ref_alignment, tree_newick, tree_sequence_alignment, tree_json
from build_tree import ancestral_sequence_inference
from Bio import Phylo, AlignIO

def tree_alignment(T):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio import AlignIO

    fname = tree_sequence_alignment(path, 'nuc')
    seqs = [SeqRecord(Seq(''.join(T.root.sequence)), name='root', id='root')]
    for node in T.find_clades():
        seqs.append(SeqRecord(Seq(''.join(node.sequence)), description='', name=node.name, id=node.name))
    return MultipleSeqAlignment(seqs)

if __name__ == '__main__':
    parser = generic_argparse("Infer mutations for one segment on the tree of another")

    args = parser.parse_args()
    path = args.path
    base_path = '/'.join(path.rstrip('/').split('/')[:-1]) + '/'

    segments = map(lambda x:x.split('/')[-1], filter(lambda x:os.path.isdir(x) and 'meta' not in x, glob.glob(base_path+'/*')))
    for segment in segments:
        T = Phylo.read(tree_newick(path), 'newick')

        sequences = base_path+ segment+'/results/ref_nuc_aln.fasta'
        try:
            aln = AlignIO.read(sequences, 'fasta')
        except:
            continue
        reference = base_path+'metadata/h3n2_%s_outgroup.gb'%segment

        tt = ancestral_sequence_inference(tree=T, aln=aln, optimize_branch_length=False, infer_gtr=True)

        selected_features = load_features(reference, feature_names=None, feature_key='gene')
        tree_aln = tree_alignment(tt.tree)

        for fname, feature in selected_features.items():
            translation = translate_feature(tree_aln, feature)
            AlignIO.write(translation, path+'/results/tree_%s_aln.fasta'%fname, 'fasta')

    T = Phylo.read(tree_newick(path), 'newick')
    tree_meta = read_tree_meta_data(path)
    for gene, aln_fname in get_genes_and_alignments(path, tree=True):
        if gene!='nuc':
            muts = get_amino_acid_mutations(T, aln_fname)
            for node_name in tree_meta:
                tree_meta[node_name][gene+'_mutations'] = muts[node_name]
    write_tree_meta_data(path, tree_meta)

    clades = {
            "3c3.a": [('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'S')],
            "3c3":   [('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'F')],
            "3c2.a": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),
                ('HA1',311,'H'), ('ha_nuc',1491,'A'), ('ha_nuc', 234, 'A')],
            "3c2.a1":[('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),
                ('HA1',311,'H'), ('ha_nuc',1491,'A'), ('HA1',171,'K'),
                ('HA2',77,'V'), ('HA2',155,'E'), ('HA2',160,'N')],
            "3c2":   [('HA1',144,'N'), ('HA1',159,'F'), ('HA1',225,'N'),
                ('HA1',160,'T'), ('HA1',142,'R')],
            "3c3.b": [('HA1',83,'R'), ('HA1',261,'Q'), ('HA1',62,'K'),
                ('HA1',122,'D')],
            "1": [('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'),
                ('ha_nuc',1491,'A'), ('ha_nuc',234,'A'), ('HA1',53,'N'),
                ('HA1',144,'R'), ('HA1',171,'K'), ('HA1',192,'T'),
                ('HA1',197,'H')],
            "2": [('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'),
                ('ha_nuc',1491,'A'), ('ha_nuc',234,'A'), ('HA1',121,'K'),
                ('HA1',144,'K')],
            "3": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),
                ('HA1',311,'H'), ('ha_nuc',1491,'A'), ('ha_nuc', 234, 'A'),
                ('HA1',131,'K'), ('HA1',142,'K'), ('HA1',261,'Q')],
            "4": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),
                ('HA1',311,'H'), ('ha_nuc',1491,'A'), ('ha_nuc', 234, 'G'),
                ('HA2',150,'E'), ('ha_nuc',114,'T')],
            "5": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'),
                ('ha_nuc',1491,'A'), ('ha_nuc', 234, 'G'), ('HA1',92,'R'),
                ('HA1',311,'Q'), ('ha_nuc',538,'C')]
        }
    tree_meta = assign_clades(path, base_path+'ha/results/tree_nuc_aln.fasta', clades)
