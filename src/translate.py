import numpy as np
import os
from filenames import ref_alignment, tree_sequence_alignment
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from util import safe_translate, load_features, generic_argparse

def translate_feature(aln, feature):
    translations = []
    for seq in aln:
        aa_seq = safe_translate(str(feature.extract(seq).seq))
        aa_seq = SeqRecord(seq=Seq(aa_seq), name=seq.name, id=seq.name)
        translations.append(aa_seq)
    return MultipleSeqAlignment(translations)


def translate(aln_fname, reference, feature_names, name_func):
    try:
        aln = AlignIO.read(aln_fname, 'fasta')
    except:
        print("Loading input alignment failed!:", aln_fname)

    selected_features = load_features(reference, feature_names)

    for fname, feature in selected_features.items():
        translation = translate_feature(aln, feature)
        AlignIO.write(translation, name_func(fname), 'fasta')


if __name__ == '__main__':
    parser = generic_argparse("Translate the nucleotide alignments")
    parser.add_argument('--reference', required=True,
                        help='genbank file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to translate")
    args = parser.parse_args()

    path = args.path

    if not args.genes:
        genes = load_features(args.reference).keys()
    else:
        genes = args.genes

    for func in [tree_sequence_alignment, ref_alignment]:
        aln_fname = func(path, 'nuc')
        if os.path.isfile(aln_fname):
            translate(aln_fname,
                  args.reference, genes,
                  lambda x:func(path, x))

