import numpy as np
from filenames import tree_newick, mugration_model
from util import read_sequence_meta_data, read_tree_meta_data, write_tree_meta_data
from util import collect_tree_meta_data, generic_argparse, TINY


def mugration_inference(tree=None, seq_meta=None, field='country', confidence=True,
                        infer_gtr=True, root_state=None, missing='?'):
        from treetime import GTR
        from Bio.Align import MultipleSeqAlignment
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq


        # Determine alphabet
        places = set()
        for meta in seq_meta.values():
            if field in meta:
                places.add(meta[field])
        if root_state is not None:
            places.add(root_state)

        # construct GTR (flat for now). The missing DATA symbol is a '-' (ord('-')==45)
        places = sorted(places)
        nc = len(places)
        if nc>180:
            print("geo_inference: can't have more than 180 places!")
            return None
        elif nc==1:
            print("geo_inference: only one place found -- set every internal node to %s!"%places[0])
            return None
        elif nc==0:
            print("geo_inference: list of places is empty!")
            return None
        else:
            alphabet = {chr(65+i):place for i,place in enumerate(places)}
            myGeoGTR = GTR.custom(pi = np.ones(nc, dtype=float)/nc, W=np.ones((nc,nc)),
                              alphabet = np.array(sorted(alphabet.keys())))
            missing_char = chr(65+nc)
            alphabet[missing_char]=missing
            myGeoGTR.profile_map[missing_char] = np.ones(nc)
            alphabet_rev = {v:k for k,v in alphabet.iteritems()}

            pseudo_seqs = []
            for name, meta in seq_meta.items():
                s=alphabet_rev[meta[field]] if field in meta else missing_char
                pseudo_seqs.append(SeqRecord(Seq(s), name=name, id=name))
            aln = MultipleSeqAlignment(pseudo_seqs)

            from treetime import TreeAnc
            tt = TreeAnc(tree=tree, aln=aln, gtr=myGeoGTR)
            tt.use_mutation_length=False
            tt.infer_ancestral_sequences(infer_gtr=infer_gtr, store_compressed=False, pc=5.0,
                                         marginal=True, normalized_rate=False)

            for node in tt.tree.find_clades():
                node.__setattr__(field, alphabet[node.sequence[0]])

            if confidence:
                for node in tt.tree.find_clades():
                    pdis = node.marginal_profile[0]
                    S = -np.sum(pdis*np.log(pdis+TINY))

                    marginal = [(alphabet[tt.gtr.alphabet[i]], pdis[i]) for i in range(len(tt.gtr.alphabet))]
                    marginal.sort(key=lambda x: x[1], reverse=True) # sort on likelihoods
                    marginal = [(a, b) for a, b in marginal if b > 0.01][:4] #only take stuff over 1% and the top 4 elements
                    conf = {a:b for a,b in marginal}
                    node.__setattr__(field + "_entropy", S)
                    node.__setattr__(field + "_confidence", conf)

            return tt, alphabet



if __name__ == '__main__':
    parser = generic_argparse("Infer ancestral states for a discrete character")
    parser.add_argument('--field', default='region',
                        help='meta data field to perform discrete reconstruction on')
    parser.add_argument('--confidence',action="store_true",
                        help='record the distribution of subleading mugration states')
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")

    args = parser.parse_args()
    path = args.path
    T = tree_newick(path)

    import time
    start = time.time()

    seq_meta = read_sequence_meta_data(path)
    tree_meta = read_tree_meta_data(path)
    end = time.time()
    print "Reading in meta files took {}".format(str(end-start))

    start = time.time()
    tt, alphabet = mugration_inference(tree=T, seq_meta=seq_meta,
                        field=args.field, confidence=args.confidence)

    end = time.time()
    print "Mugration took {}".format(str(end-start))

    if args.confidence:
        fields = [args.field, args.field+'_confidence', args.field + '_entropy']
    else:
        fields = [args.field]

    start = time.time()
    collect_tree_meta_data(tt.tree, fields, isvcf=args.vcf, meta=tree_meta)
    write_tree_meta_data(path, tree_meta)

    with open(mugration_model(path, args.field),'w') as ofile:
        ofile.write('Map from character to field name\n')
        for k,v in alphabet.items():
            ofile.write(k+':\t'+v+'\n')
        ofile.write('\n\n')

        ofile.write(str(tt.gtr))

    end = time.time()
    print "Mugration file writing took {}".format(str(end-start))