import numpy as np
from Bio import AlignIO, Seq
from filenames import raw_alignment, ref_alignment
from util import write_fasta, generic_argparse

def strip_non_reference(path, reference, keep_reference=False):
    '''
    return sequences that have all insertions relative to the reference
    removed. The alignment is read from file and returned as list of sequences.
    '''
    aln = AlignIO.read(raw_alignment(path), 'fasta')
    seqs = {s.name:s for s in aln}
    if reference in seqs:
        ref_array = np.array(seqs[reference])
        ungapped = ref_array!='-'
        ref_aln_array = np.array(aln)[:,ungapped]
    else:
        print("reference", reference, "not found in alignment")
        return

    out_seqs = []
    for seq, seq_array in zip(aln, ref_aln_array):
        seq.seq = Seq.Seq(''.join(seq_array))
        if keep_reference or seq.name!=reference:
            out_seqs.append(seq)

    return out_seqs

if __name__ == '__main__':
    parser = generic_argparse("strip out all positions that don't align to the reference")
    parser.add_argument('--reference', required=True,
                        help='the name of the reference sequence')
    parser.add_argument('--keep_reference', action='store_true', default=False,
                        help='keep the reference as part of the alignment')
    args = parser.parse_args()

    seqs = strip_non_reference(args.path, args.reference)
    write_fasta(seqs, ref_alignment(args.path))
