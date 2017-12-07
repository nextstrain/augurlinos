import pandas as pd
import numpy as np
import os, shutil
from filenames import ref_alignment, tree_newick, sequence_gtr_model, tree_sequence_alignment
from util import read_sequence_meta_data, parse_date, write_tree_meta_data
from util import collect_tree_meta_data, write_json, generic_argparse
from Bio import Phylo


def build_fasttree(aln_file, out_file, clean_up=True):
    call = ["fasttree", "-nt", aln_file, "1>", out_file, "2>", "fasttree.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    try:
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove('fasttree.log')
    except:
        print("TREE BUILDING FAILED")
        T=None

    return T


def build_iqtree(aln_file, out_file, clean_up=True, nthreads=2):
    call = ["iqtree", "-nt", nthreads, "-s", aln_file, ">", "iqtree.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    try:
        T = Phylo.read(aln_file+".treefile", 'newick')
        if clean_up:
            os.remove('iqtree.log')
            os.remove(aln_file+".*")
    except:
        print("TREE BUILDING FAILED")
        T=None
    return T


def timetree(tree=None, aln=None, seq_meta=None, keeproot=False,
             confidence=False, resolve_polytomies=True, max_iter=2,
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, **kwarks):
    from treetime import TreeTime
    dates = {}
    for name, data in seq_meta.items():
        num_date = parse_date(data["date"], date_fmt)
        if num_date is not None:
            dates[name] = num_date

    tt = TreeTime(tree=tree, aln=aln, dates=dates, gtr='JC69')

    if confidence and use_marginal:
        # estimate confidence intervals via marginal ML and assign marginal ML times to nodes
        marginal = 'assign'
    else:
        marginal = confidence

    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
           resolve_polytomies=resolve_polytomies, max_iter=max_iter, **kwarks)

    for n in T.find_clades():
        n.num_date = n.numdate # treetime convention is different from augur...
        # get 90% max posterior region)
        if confidence:
            n.num_date_confidence = list(tt.get_max_posterior_region(n, 0.9))
    return tt


def ancestral_sequence_inference(tree=None, aln=None, infer_gtr=True,
                                 optimize_branch_length=True):
    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, gtr='JC69')

    if optimize_branch_length:
        tt.optimize_seq_and_branch_len(infer_gtr=infer_gtr)
    else: # only infer ancestral sequences, leave branch length untouched
        tt.infer_ancestral_sequences(infer_gtr=infer_gtr)

    return tt


def export_sequence_fasta(T, path):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    from Bio import AlignIO

    fname = tree_sequence_alignment(path, 'nuc')
    seqs = [SeqRecord(Seq(''.join(T.root.sequence)), name='root', id='root')]
    for node in T.find_clades():
        seqs.append(SeqRecord(Seq(''.join(node.sequence)), name=node.name, id=node.name))
    AlignIO.write(MultipleSeqAlignment(seqs), fname, 'fasta')


if __name__ == '__main__':
    parser = generic_argparse("Build the tree from the prepared sequence data")
    parser.add_argument('--nthreads', type=int, default=2,
                        help='number of threads')
    parser.add_argument('--ancestral', action='store_true', default=False,
                        help='calculate ancestral sequences')
    parser.add_argument('--timetree', action='store_true', default=False,
                       help='infer time stamped phylogeny')
    parser.add_argument('--confidence', action='store_true', default=False,
                       help='estimate confidence intervals for node timing')
    parser.add_argument('--Tc', type=float, default=0.0,
                       help='coalescence time scale measured in substitution rate units')
    parser.add_argument('--keeproot', action='store_true', default=False,
                        help="don't reroot the tree")
    args = parser.parse_args()
    path = args.path

    date_fmt = '%Y-%m-%d'

    T = build_fasttree(ref_alignment(path), tree_newick(path))
    #T = tree_newick(path)
    meta = read_sequence_meta_data(path)
    fields = ['branchlength', 'clade']

    if args.timetree:
        tt = timetree(tree=T, aln=ref_alignment(path), confidence=args.confidence,
                      seq_meta=meta, reroot=None if args.keeproot else 'best',Tc=args.Tc)
        T = tt.tree
        fields.extend(['mutations', 'mutation_length', 'num_date', 'clock_length'])
        if args.confidence:
            fields.append('num_date_confidence')
    elif args.ancestral:
        tt = ancestral_sequence_inference(tree=T, aln=ref_alignment(path))
        T = tt.tree
        fields.extend(['mutations', 'mutation_length'])

    clade_index = 0
    for n in T.find_clades(order='preorder'):
        n.clade = clade_index
        clade_index+=1

    Phylo.write(T, tree_newick(path), 'newick')
    meta_dic = collect_tree_meta_data(T, fields)
    write_tree_meta_data(path, meta_dic)

    with open(sequence_gtr_model(path),'w') as ofile:
        ofile.write(str(tt.gtr))

    if args.timetree or args.ancestral:
        export_sequence_fasta(T, path)
