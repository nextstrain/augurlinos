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


def build_iqtree(aln_file, out_file, iqmodel, clean_up=True, nthreads=2):
    #return Phylo.read(out_file.replace(".nwk",".iqtree.nwk"), 'newick') #uncomment for debug skip straight to TreeTime
    if iqmodel:
        call = ["iqtree", "-nt", str(nthreads), "-s", aln_file, "-m", iqmodel[0],
            ">", "iqtree.log"]
    else:
        call = ["iqtree", "-nt", str(nthreads), "-s", aln_file, ">", "iqtree.log"]
    print(" ".join(call))
    os.system(" ".join(call))
    try:
        T = Phylo.read(aln_file+".treefile", 'newick')
        shutil.copyfile(aln_file+".treefile", out_file)
        #this allows the user to check intermediate output, as tree.nwk will be
        #written over with TreeTime tree
        shutil.copyfile(aln_file+".treefile", out_file.replace(".nwk",".iqtree.nwk"))
        if clean_up:
            #allow user to see chosen model
            shutil.copyfile('iqtree.log', out_file.replace("tree.nwk","iqtree.log"))
            os.remove('iqtree.log')
            for filename in glob.glob(aln_file+".*"):
                os.remove(filename)
    except:
        print("TREE BUILDING FAILED")
        T=None
    return T


def timetree(tree=None, aln=None, ref=None, seq_meta=None, keeproot=False,
             confidence=False, resolve_polytomies=True, max_iter=2, dateLimits=None,
             infer_gtr=True, Tc=0.01, reroot='best', use_marginal=False, **kwarks):
    from treetime import TreeTime

    dL_int = None
    if dateLimits:
        dL_int = [int(x) for x in dateLimits]
        dL_int.sort()

    dates = {}
    for name, data in seq_meta.items():
        num_date = parse_date(data["date"], date_fmt, dL_int)
        if num_date is not None:
            dates[name] = num_date

    #send ref, if is None, does no harm
    tt = TreeTime(tree=tree, aln=aln, ref=ref, dates=dates, gtr='JC69')

    if confidence and use_marginal:
        # estimate confidence intervals via marginal ML and assign marginal ML times to nodes
        marginal = 'assign'
    else:
        marginal = confidence

    #Length of VCF files means GTR model with gaps causes overestimation of mutation TO gaps
    #so gaps appear in internal nodes when no gaps at tips! To prevent....
    pi = None
    if ref != None: #if VCF, fix pi
        pi = np.array([0.1618, 0.3188, 0.3176, 0.1618, 0.04]) #from real runs (Walker 2018)


    tt.run(infer_gtr=infer_gtr, root=reroot, Tc=Tc, time_marginal=marginal,
        resolve_polytomies=resolve_polytomies, max_iter=max_iter, fixed_pi=pi, **kwarks)





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

def export_sequence_VCF(tt, path, var_ambigs=False):
    tree_dict = tt.get_tree_dict(var_ambigs)
    #in new augur, set 'compress' based on input file ending!
    write_VCF_style_alignment(tree_dict, tree_vcf_alignment(path,'nuc'), compress=False)


def write_out_variable_fasta(compress_seq, path, drmfile=None):
    #Now also writes out a file tracking the actual position of the variable site written out.
    #So if viewing the variable Fasta and wondering what real position column 5 is,
    #one can look at row 5 of this file, and it will give the real bp position.

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    sequences = compress_seq['sequences']
    ref = compress_seq['reference']
    positions = compress_seq['positions']

    #If want to exclude DRMs from initial treebuild, read in here
    #and remove from the positions we're going to include in var sites
    if drmfile:
        from util import read_in_DRMs
        DRM_info = read_in_DRMs(drmfile)
        drmPositions = DRM_info['drmPositions']

        toDel = []
        for drmKey in drmPositions:
            temp = np.where(positions==drmKey)[0]
            if len(temp)!=0:
                toDel.append(temp[0])

        positions = np.delete(positions, toDel)

    #get sequence names
    seqNames = sequences.keys()

    #get the variables sites, either ALT or REF as already determined above
    #use faster method
    sites = []
    pos = []
    for key in positions:
        pattern = []
        for k,v in sequences.iteritems():
            try:
                pattern.append(sequences[k][key])
            except KeyError, e:
                pattern.append(ref[key])
        #pattern = [ sequences[k][key] if key in sequences[k].keys() else ref[key] for k,v in sequences.iteritems() ]
        origPattern = list(pattern)
        if '-' in pattern:
            pattern = [value for value in origPattern if value != '-']
            #remove gaps to see if otherwise non-variable
            #print pattern
        un = np.unique(pattern, return_counts=True)
        if len(un[0])==0 or len(un[0])==1 or (len(un[0])==2 and min(un[1])==1): #'singleton' mutation - only happens in 1 seq
        #if len(un[0])==1: #don't write out identical bases
            False #don't append! (python makes me put something here)
        else:
            sites.append(origPattern) #append original with gaps (if present)
            pos.append(str(key))

    #rotate into an alignment and turn into list of SeqRecord to output easily
    sites = np.asarray(sites)
    align = np.rot90(sites)
    seqNamesCorr = list(reversed(seqNames))
    toFasta = [ SeqRecord(id=seqNamesCorr[i], seq=Seq("".join(align[i])), description='') for i in xrange(len(sequences.keys()))]

    with open(var_site_alignment(path)+".positions.txt", 'w') as the_file:
        the_file.write("\n".join(pos))

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
    parser.add_argument('--fasttree', action='store_true', default=True,
                        help="use fasttree for initial tree building (default)")
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")

    #EBH 5 Jan 2018
    parser.add_argument('--iqmodel', nargs=1, help='model to use with iqtree')
    parser.add_argument('--drm', type=str,
                        help="file of DRMs to exclude from inital tree-building")

    #EBH 14 Feb 2018
    parser.add_argument('--roottype', nargs="+",#type=str, default="residual",
                        help="type of rerooting. options are 'rsq', 'residual' (default), and 'oldest'")

    #EBH 16 Mar 2018
    parser.add_argument('--varAmbigs', action='store_true', default=False,
                        help="preserve ambiguous bases at variable sites in recorded mutations")
    parser.add_argument('--dateLimit', nargs='+',
                        help="specify min and max year for samples without dates. Order doesn't matter. If only one value, taken as min, and max set to current year.")

    args = parser.parse_args()
    path = args.path

    date_fmt = '%Y-%m-%d'

    import time

    if args.vcf:
        #read in VCF and reference fasta and store
        start = time.time()
        compress_seq = read_in_vcf(recode_gzvcf_name(path), ref_fasta(path))
        sequences = compress_seq['sequences']
        ref = compress_seq['reference']
        end = time.time()
        print "Reading in VCF took {}".format(str(end-start))

        #write out the reduced fasta (only variable sites) to be used by
        #by iqtree/raxml/fasttree  ("var_site_alignment(path)")
        start = time.time()
        write_out_variable_fasta(compress_seq, path, args.drm)
        end = time.time()
        print "Writing out variable sites took {}".format(str(end-start))

    if args.iqmodel and not args.iqtree:
        print "Cannot specify model unless using IQTree. Model specification ignored."

    if args.vcf:
        treebuild_align = var_site_alignment(path)
    else:
        treebuild_align = ref_alignment(path)

    start = time.time()
    if args.raxml:
        T = build_raxml(treebuild_align, tree_newick(path), path, args.nthreads)
    elif args.iqtree:
        T = build_iqtree(treebuild_align, tree_newick(path), args.iqmodel, args.nthreads)
    else: #use fasttree - if add more options, put another check here
        T = build_fasttree(treebuild_align, tree_newick(path))
    end = time.time()
    print "Building original tree took {}".format(str(end-start))

    meta = read_sequence_meta_data(path)
    fields = ['branchlength', 'clade']

    #Anything but a list of sequences to root by, shouldn't go as a "list".
    if len(args.roottype) == 1:
        args.roottype = args.roottype[0]

    start = time.time()
    if args.timetree:
        if args.vcf:
            tt = timetree(tree=T, aln=sequences, ref=ref, confidence=args.confidence, dateLimits=args.dateLimit,
                          seq_meta=meta, reroot=None if args.keeproot else args.roottype, Tc=args.Tc)#, use_marginal=True)
        else:
            tt = timetree(tree=T, aln=ref_alignment(path), confidence=args.confidence, dateLimits=args.dateLimit,
                          seq_meta=meta, reroot=None if args.keeproot else args.roottype, Tc=args.Tc)

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

    end = time.time()
    print "TreeTime took {}".format(str(end-start))

    clade_index = 0
    for n in T.find_clades(order='preorder'):
        n.clade = clade_index
        clade_index+=1

    if args.vcf:
        Phylo.write(T, tree_newick(path), 'newick', format_branch_length='%1.8f')
        #This gives more digits to branch length which for small numbers
        #means it's more than just 0 and 0.00001 in the Newick!
    else:
        Phylo.write(T, tree_newick(path), 'newick')

    if args.varAmbigs: #if requested, put ambigs back on tips
        tt.recover_var_ambigs()
    meta_dic = collect_tree_meta_data(T, fields, args.vcf)
    write_tree_meta_data(path, meta_dic)

    with open(sequence_gtr_model(path),'w') as ofile:
        ofile.write(str(tt.gtr))

    start = time.time()
    #do NOT print out all full sequences if VCF - will be huge!
    if args.timetree or args.ancestral:
        if args.vcf:
            export_sequence_VCF(tt, path, args.varAmbigs)
        else:
            export_sequence_fasta(T, path)
    end = time.time()
    print "Writing out VCF/Fasta took {}".format(str(end-start))