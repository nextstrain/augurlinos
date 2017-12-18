import numpy as np
import os
from filenames import (ref_alignment, tree_sequence_alignment, tree_vcf_alignment,
                        recode_gzvcf_name, dropped_genes, translation_vcf_file, translation_ref_file,
                        ref_fasta, tree_newick)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from util import safe_translate, load_features, generic_argparse, read_in_vcf, write_VCF_translation
from util import read_tree_meta_data, write_tree_meta_data, get_genes_and_alignments

def translate_feature(aln, feature):
    translations = []
    for seq in aln:
        aa_seq = safe_translate(str(feature.extract(seq).seq))
        aa_seq = SeqRecord(seq=Seq(aa_seq), name=seq.name, id=seq.name)
        translations.append(aa_seq)
    return MultipleSeqAlignment(translations)


def translate(aln_fname, reference, name_func, feature_names=None):
    try:
        aln = AlignIO.read(aln_fname, 'fasta')
    except:
        print("Loading input alignment failed!:", aln_fname)

    selected_features = load_features(reference, feature_names)

    for fname, feature in selected_features.items():
        translation = translate_feature(aln, feature)
        AlignIO.write(translation, name_func(fname), 'fasta')


def translate_vcf_feature(sequences, ref, feature):
    prot = {}
    prot['sequences'] = {}
    prot['positions'] = []

    refNuc = str(feature.extract( SeqRecord(seq=Seq(ref)) ).seq)
    ref_aa_seq = safe_translate(refNuc)
    prot['reference'] = ref_aa_seq

    start = int(feature.location.start)
    end = int(feature.location.end)

    for seqk in sequences.keys():
        #get positions where diffs
        varSite = np.array(sequences[seqk].keys())
        #reduce to only those within current gene
        geneVarSites = np.logical_and(varSite >= start, varSite <= end)
        #translate this back to nuc position
        nucVarSites = varSite[geneVarSites]
        #get it in position within the gene! - because whole genome may not be in frame!! But we must assume gene is..
        genNucSites = nucVarSites-start

        #Translate the codon this nuc diff is in, and find out which AA loc
        #But need numbering to be w/in protin, not whole genome!
        aaRepLocs = {i//3:safe_translate( "".join([sequences[seqk][key+start]
                                if key+start in sequences[seqk].keys() else refNuc[key]
                            for key in range(i-i%3,i+3-i%3)]) )
                        for i in genNucSites}

        aaRepLocsFinal = {}
        #remove if is a synonymous mutation!
        for key,val in aaRepLocs.iteritems():
            if ref_aa_seq[key] != val:
                aaRepLocsFinal[key] = val
                #print "new codon ({}) is same as ref({}) at position {}".format(val, ref_aa_seq_np[key], key)

        aaRepLocs = aaRepLocsFinal

        #store the dict of differences
        prot['sequences'][seqk] = aaRepLocs

        #add to list of positions if needed
        for key in aaRepLocs.keys():
            if key not in prot['positions']:
                prot['positions'].append(key)

        #IF we wanted to get full translation....
        #ref_aa_seq_np = np.array(list(ref_aa_seq), 'S1')
        #newAA = ref_aa_seq_np.copy()
        #newAA[aaRepLocs.keys()] = aaRepLocs.values()

    prot['positions'].sort()

    #if no variable sites, exclude this gene
    if len(prot['positions']) == 0:
        return None
    else:
        return prot


def translate_vcf(vcf_fname, reference, path, feature_names=None):
    try:
        vcf_dict = read_in_vcf(vcf_fname, ref_fasta(path), compressed=False )
    except:
        print "Loading input alignment failed!: {}".format(vcf_fname)

    selected_features = load_features(reference, feature_names)

    ref = vcf_dict['reference']
    sequences = vcf_dict['sequences']

    prots = {}
    deleted = []
    notMult3 = []

    import time
    start = time.time()
    #if genes have no mutations across sequences, they are dropped here from further analysis
    #check that gene lengths are multiples of 3. The first occurance causes an error in
    #Biopython, but subsequent ones do not - so make it ourselves.
    for fname,feature in selected_features.items():
        if len(str(feature.extract( SeqRecord(seq=Seq(ref)) ).seq))%3 != 0:
            notMult3.append(fname)

        prot_dict = translate_vcf_feature(sequences, ref, feature)
        if prot_dict is not None:
            prots[fname] = prot_dict
        else:
            deleted.append(fname)
    end = time.time()
    print "Translations took {}".format(str(end-start))

    start = time.time()
    #print out VCF of proteins
    write_VCF_translation(prots, translation_vcf_file(path), translation_ref_file(path))
    end = time.time()
    print "Writing out VCF took {}".format(str(end-start))

    #warn of those that don't have a length mult of 3
    print "WARNING: These genes do not have lengths that are a multiple of 3!\n{}".format(str(notMult3))

    #print dropped genes to a text file
    if len(deleted) != 0:
        with open(dropped_genes(path), 'w') as the_file:
            for d in deleted:
                the_file.write(d+"\n")
        print "{} genes had no mutations and so will be excluded. Excluded genes can be found in {}".format(len(deleted), dropped_genes(path))

    return prots


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


def assign_amino_acid_muts_vcf(prots, path):
    tree_meta = read_tree_meta_data(path)
    seqNames = prots[prots.keys()[0]]['sequences'].keys()

    #go through every gene in the prots nested dict
    for fname, prot in prots.iteritems():
        sequences = prot['sequences']
        ref = prot['reference']
        positions = prot['positions']

        pats = []
        i=0
        #for each position, get the mutation in the right format
        #[ancestral][position][mutation]
        while i < len(positions):
            pi = positions[i]
            refb = ref[pi]

            pattern = [ refb+str(pi)+sequences[k][pi] if pi in sequences[k].keys()
                        else "" for k,v in sequences.iteritems() ]
            pats.append(pattern)
            i+=1

        #convert our list of lists to matrix
        patMat = np.matrix(pats)

        #for every sequence, assign the mutations in tree_meta
        for i in xrange(len(seqNames)):
            node_name = seqNames[i]
            ary = np.array(patMat[:,i]).reshape(-1,)
            tree_meta[node_name][fname+'_mutations'] = ",".join(ary[ary != ''])

    #write it out!
    write_tree_meta_data(path, tree_meta)


def get_genes_from_file(fname):
    #based on get_dropped_strains from prepare.vcf
    genes = []
    if os.path.isfile(fname):
        with open(fname) as ifile:
            for line in ifile:
                fields = line.strip().split('#')
                if fields[0].strip():
                    genes.append(fields[0].strip())
    else:
        print("File with genes not found. Looking for", fname)

    return genes


if __name__ == '__main__':
    parser = generic_argparse("Translate the nucleotide alignments")
    parser.add_argument('--reference', required=True,
                        help='genbank file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to translate")
    #EBH 11 Dec 17
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")
    parser.add_argument('--assignMuts', action='store_true', default=False,
                        help="write amino acid mutations onto the tree")
    args = parser.parse_args()

    path = args.path

    #The original way of doing this called load_features twice!
    if not args.genes:
        genes = None #if load_features is passed None it loads all
    else:
        genes = args.genes
        if os.path.isfile(genes[0]):
            genes = get_genes_from_file(genes[0])

    if args.vcf:
        #VCF has no 'ref_' alignment as this doesn't make sense
        #use the tree alignment if there, otherwise use the original VCF
        inputAlign = recode_gzvcf_name(path)
        if os.path.isfile(tree_vcf_alignment(path)):
            inputAlign = tree_vcf_alignment(path)

        prots = translate_vcf(inputAlign, args.reference, path, genes)

        if args.assignMuts:
            assign_amino_acid_muts_vcf(prots, path)

    else:
        for func in [tree_sequence_alignment, ref_alignment]:
            aln_fname = func(path, 'nuc')
            if os.path.isfile(aln_fname):
                translate(aln_fname,
                      args.reference,
                      lambda x:func(path, x), genes)

        if args.assignMuts:
            tree_meta = read_tree_meta_data(path)
            T = Phylo.read(tree_newick(path), 'newick')

            for gene, aln_fname in get_genes_and_alignments(path, tree=True):
                if gene!='nuc':
                    muts = get_amino_acid_mutations(T, aln_fname)

                for node_name in tree_meta:
                    tree_meta[node_name][gene+'_mutations'] = muts[node_name]
            write_tree_meta_data(path, tree_meta)