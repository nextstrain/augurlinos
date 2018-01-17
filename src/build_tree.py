import pandas as pd
import numpy as np
import os, shutil, glob
from filenames import ref_alignment, tree_newick, sequence_gtr_model, tree_sequence_alignment
from filenames import recode_gzvcf_name, ref_fasta, var_site_alignment, tree_vcf_alignment
from util import read_sequence_meta_data, parse_date, write_tree_meta_data, write_VCF_style_alignment
from util import collect_tree_meta_data, write_json, generic_argparse, read_in_vcf
from Bio import Phylo


def build_raxml(aln_file, out_file, path, clean_up=True, nthreads=2):
    call = ["raxml","-f d -T",str(nthreads),"-m GTRCAT -c 25 -p 235813 -n tre -s",aln_file,"-w",os.getcwd()+"/"+path+"/results/","> raxml.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    shutil.copy(path+'/results/RAxML_bestTree.tre', out_file)
    try:
        T = Phylo.read(out_file, 'newick')
        if clean_up:
            os.remove('raxml.log')
    except:
        print("TREE BUILDING FAILED")
        T=None

    return T

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
    call = ["iqtree", "-nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    try:
        T = Phylo.read(aln_file+".treefile", 'newick')
        shutil.copyfile(aln_file+".treefile", out_file)
        if clean_up:
            os.remove('iqtree.log')
            for filename in glob.glob(aln_file+".*"):
                os.remove(filename)
    except:
        print("TREE BUILDING FAILED")
        T=None
    return T


def timetree(tree=None, aln=None, ref=None, seq_meta=None, keeproot=False,
             confidence=False, resolve_polytomies=True, max_iter=2,
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, **kwarks):
    from treetime import TreeTime
    dates = {}
    for name, data in seq_meta.items():
        num_date = parse_date(data["date"], date_fmt)
        if num_date is not None:
            dates[name] = num_date

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates, gtr='JC69')

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


def ancestral_sequence_inference(tree=None, aln=None, ref=None, infer_gtr=True,
                                 optimize_branch_length=True):
    from treetime import TreeAnc
    tt = TreeAnc(tree=tree, aln=aln, ref=ref, gtr='JC69')

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
        seqs.append(SeqRecord(Seq(''.join(node.sequence)), description='', name=node.name, id=node.name))
    AlignIO.write(MultipleSeqAlignment(seqs), fname, 'fasta')

def export_sequence_VCF(tt, path):
    tree_dict = tt.get_tree_dict()
    write_VCF_style_alignment(tree_dict, tree_vcf_alignment(path,'nuc'))


def write_out_variable_fasta(compress_seq, path):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    sequences = compress_seq['sequences']
    ref = compress_seq['reference']
    positions = compress_seq['positions']

    #get sequence names
    seqNames = sequences.keys()

    #get the variables sites, either ALT or REF as already determined above
    sites = []
    for key in positions:
        pattern = [ sequences[k][key] if key in sequences[k].keys() else ref[key] for k,v in sequences.iteritems() ]
        sites.append(pattern)

    #rotate into an alignment and turn into list of SeqRecord to output easily
    sites = np.asarray(sites)
    align = np.rot90(sites)
    toFasta = [ SeqRecord(id=seqNames[i], seq=Seq("".join(align[i])), description='') for i in xrange(len(sequences.keys()))]

    #now output this as fasta to read into raxml or iqtree
    SeqIO.write(toFasta, var_site_alignment(path), 'fasta')


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

    #EBH 4 Dec 2017
    parser.add_argument('--iqtree', action='store_true', default=False,
                        help="use iqtree for initial tree building")
    parser.add_argument('--raxml', action='store_true', default=False,
                        help="use raxml for initial tree building")
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")

    args = parser.parse_args()
    path = args.path

    date_fmt = '%Y-%m-%d'

    if args.vcf:
        #read in VCF and reference fasta and store
        compress_seq = read_in_vcf(recode_gzvcf_name(path), ref_fasta(path))
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']

        #write out the reduced fasta (only variable sites) to be read in
        #by iqtree/raxml/fasttree  ("var_site_alignment(path)")
        write_out_variable_fasta(compress_seq, path)

    if args.vcf:
        treebuild_align = var_site_alignment(path)
    else:
        treebuild_align = ref_alignment(path)

    if args.raxml:
        T = build_raxml(treebuild_align, tree_newick(path), path)
    elif args.iqtree:
        T = build_iqtree(treebuild_align, tree_newick(path))
    else: #use fasttree
        T = build_fasttree(treebuild_align, tree_newick(path), clean_up=False)

    meta = read_sequence_meta_data(path)
    fields = ['branchlength', 'clade']

    if args.timetree:
        if args.vcf:
            tt = timetree(tree=T, aln=sequences, ref=ref, confidence=args.confidence,
                          seq_meta=meta, reroot=None if args.keeproot else 'best', Tc=args.Tc)
        else:
            tt = timetree(tree=T, aln=ref_alignment(path), confidence=args.confidence,
                          seq_meta=meta, reroot=None if args.keeproot else 'best', Tc=args.Tc)

        T = tt.tree
        fields.extend(['mutations', 'mutation_length', 'num_date', 'clock_length'])
        if args.confidence:
            fields.append('num_date_confidence')
    elif args.ancestral:
        if args.vcf:
            tt = ancestral_sequence_inference(tree=T, aln=sequences, ref=ref)
        else:
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

    #do NOT print out all full sequences if VCF - will be huge!
    if args.timetree or args.ancestral:
        if args.vcf:
            export_sequence_VCF(tt, path)
        else:
            export_sequence_fasta(T, path)
