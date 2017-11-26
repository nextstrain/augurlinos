from frequency_estimators import alignment_frequencies, make_pivots
from filenames import ref_sequence_alignment
from util import (read_meta_data, get_genes_and_alignments, generic_argparse)


def filter_alignment(aln, meta, region=None, lower_tp=None, upper_tp=None):
    from Bio.Align import MultipleSeqAlignment
    tmp = aln
    if region is not None:
        if type(region)==str:
            tmp = [s for s in tmp if meta[s.name]['region']==region]
        elif type(region)==list:
            tmp = [s for s in tmp if meta[s.name]['region'] in region]
        else:
            self.log.warn("region must be string or list")
            return
    if lower_tp is not None:
        tmp = [s for s in tmp if np.mean(meta[s.name]['numdate'])>=lower_tp]
    if upper_tp is not None:
        tmp = [s for s in tmp if np.mean(meta[s.name]['numdate'])<upper_tp]
    return MultipleSeqAlignment(tmp)


def estimate_mutation_frequencies(aln, tps, seq_type='aa',
                                  inertia=0.0,
                                  min_freq=0.01,
                                  stiffness=20.0,
                                  in_pivots=24,
                                  include_set={}):
    '''
    calculate the frequencies of mutation
    '''
    pivots = make_pivots(in_pivots, tps)

    # instantiate alignment frequency
    aln_frequencies = alignment_frequencies(tmp_aln, tps, pivots,
                                    ws=max(2,len(tps)//10),
                                    inertia=inertia,
                                    stiffness=stiffness, method='SLSQP')
    if seq_type=='nuc': # if this is a nucleotide alignment, set all non-canonical states to N
        A = aln_frequencies.aln
        A[~((A=='A')|(A=='C')|(A=='G')|(A=='T')|('A'=='-'))] = 'N'

    aln_frequencies.mutation_frequencies(min_freq=min_freq, include_set=tmp_include_set,
                                         ignore_char='N' if seq_type=='nuc' else 'X')

    return aln_frequencies.frequencies, aln_frequencies.calc_confidence(), aln_frequencies.counts


if __name__ == '__main__':
    parser = generic_argparse("Assign amino acid mutations to the tree")
    args = parser.parse_args()
    path = args.path

    meta_dic = read_meta_data(path)
    for gene, aln_fname in get_genes_and_alignments(path, tree=False):
    	aln = AlignIO.read(aln_fname, 'fasta')
    	for region in ['global', 'NA', 'EU']:
    		tmp_aln =filter_alignment(aln, meta, region=region)

    	tps = np.array([np.mean(meta_dic[x.name]['numdate']) for x in tmp_aln])
    	freq, conf, counts = estimate_mutation_frequencies(tmp_aln, tps,
    							seq_type='nuc' if gene=='nuc' else 'aa')

