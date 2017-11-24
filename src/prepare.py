from __future__ import division, print_function
import os
from Bio import SeqIO
import pandas as pd
from filenames import meta_file_name, sequence_input
from util import write_fasta, write_sequence_meta_data, generic_argparse

def parse_fasta(fname, fields, sep='|'):
    '''
    parse a fasta file and separate the meta data from the header
    '''
    seqs = SeqIO.parse(fname, format='fasta')
    sequences = {}
    meta_data = {}
    for seq in seqs:
        tmp_meta = seq.name.split(sep)
        seq_name = tmp_meta[0]
        seq.id = seq.name = seq_name
        if seq_name in sequences:
            print("NAME CLASH:", seq_name, "found more than once")
        else:
            sequences[seq_name] = seq
            meta_data[seq_name] = {k:tmp_meta[i] for i,k in fields.items() if i<len(tmp_meta)}

    return sequences, meta_data


if __name__ == '__main__':
    parser = generic_argparse("parse fasta file and separate meta_data into table")
    parser.add_argument("--sequences", required=True, type=str,
                        help = "file with input sequences as fasta")
    args = parser.parse_args()

    if not args.path:
        path = '.'.join(os.path.basename(args.sequences).split('.')[:-1])
    else:
        path = args.path

    header_fields = {0:'strain', 2:'accession', 3:'date', 4:'region', 5:'country',
                    6:'division', 8:'db', 10:'authors', 11:'url', 12:'title',
                    13: 'journal', 14: 'paper_url'}

    sequences, meta = parse_fasta(args.sequences, header_fields)

    meta_data = pd.DataFrame(meta.values(), columns = [header_fields[i]
                                    for i in sorted(header_fields.keys())])


    write_sequence_meta_data(path, meta_data)
    write_fasta(sequences.values(), sequence_input(path), ungap=True)
