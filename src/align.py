import os
from filenames import sequence_input, raw_alignment
from util import generic_argparse

if __name__ == '__main__':
    parser = generic_argparse("Align sequences")
    parser.add_argument('--nthreads', type=int, default=2,
                        help="number of threads used by mafft")
    parser.add_argument('--aligner', default='mafft',
                        help="analysis path, e.g. zika")
    args = parser.parse_args()

    in_file = sequence_input(args.path)
    out_file = raw_alignment(args.path)

    if args.aligner=='mafft':
        os.system("mafft --anysymbol --thread %d %s 1> %s 2>mafft_stderr"%(args.nthreads, in_file, out_file))
    else:
        print('not implemented')
