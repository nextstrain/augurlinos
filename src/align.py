import os
from filenames import sequence_input, raw_alignment
from util import generic_argparse, load_reference

if __name__ == '__main__':
    parser = generic_argparse("Align sequences")
    parser.add_argument('--nthreads', type=int, default=2,
                        help="number of threads used by mafft")
    parser.add_argument('--aligner', default='mafft',
                        help="analysis path, e.g. zika")
    parser.add_argument('--reference', default=None,
                        help="reference sequence to include in alignment")
    args = parser.parse_args()

    in_file = sequence_input(args.path)
    out_file = raw_alignment(args.path)

    ref_added=False
    if args.reference and os.path.isfile(args.reference):
        from Bio import SeqIO
        ref = load_reference(args.reference)
        seqs = list(SeqIO.parse(in_file, 'fasta'))
        seq_names = [seq.id for seq in seqs]
        if ref.name not in seq_names:
            seqs.append(ref)
            ref_added = True
        else:
            ref_ii = seq_names.index(ref.id)
            seqs[ref_ii]=ref
        SeqIO.write(seqs, in_file, 'fasta')

    if args.aligner=='mafft':
        os.system("mafft --anysymbol --thread %d %s 1> %s 2>mafft_stderr"%(args.nthreads, in_file, out_file))
    else:
        print('not implemented')

    from Bio import AlignIO
    aln = AlignIO.read(out_file, 'fasta')
    for seq in aln:
        seq.seq = seq.seq.upper()

    AlignIO.write(aln, out_file, 'fasta')
