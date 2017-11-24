import numpy as np

def generic_argparse(desc):
    import argparse
    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--path', default='',
                        help="Prefix for output files (e.g., 'zika' for output like 'auspice/zika_meta.json')")
    return parser


def write_fasta(seqs, fname, ungap=False):
    from Bio import SeqIO
    with open(fname, 'w') as ofile:
        for seq in seqs:
            seq.description=""
            seq.name=seq.id
            if ungap:
                seq.seq.ungap('-')
            SeqIO.write(seq, ofile, format='fasta')


def read_sequence_meta_data(path):
    from filenames import meta_file_name
    import pandas as pd
    df = pd.read_csv(meta_file_name(path), sep='\t')
    return {m[0]:m.to_dict() for mi, m in df.iterrows()}


def write_sequence_meta_data(path, df):
    from filenames import meta_file_name
    df.to_csv(meta_file_name(path), sep='\t', index=False)


def read_tree_meta_data(path):
    from filenames import tree_meta_file_name
    import pandas as pd
    df = pd.read_csv(tree_meta_file_name(path), sep='\t')
    return {m[0]:m.to_dict() for mi, m in df.iterrows()}


def write_tree_meta_data(path, meta_dic):
    from filenames import tree_meta_file_name
    import pandas as pd
    header = set()
    for m in meta_dic.values():
        header.update(m.keys())
    header = sorted(header, key=lambda x:x!='name')
    df = pd.DataFrame(meta_dic.values(), columns = header)
    df.to_csv(tree_meta_file_name(path), sep='\t', index=False)


def collect_tree_meta_data(T, fields, meta=None):
    def mutation_format(muts):
        return ",".join(['%s%d%s'%(x[0], x[1], x[2]) for x in muts])

    if meta is None:
        meta = {}
    for n in T.find_clades():
        meta_dic = {'name':n.name}
        for field in fields:
            if hasattr(n,field):
                if 'mutations' in field:
                    meta_dic[field]=mutation_format(n.__getattribute__(field))
                else:
                    meta_dic[field]=n.__getattribute__(field)

        if n.name in meta:
            meta[n.name].update(meta_dic)
        else:
            meta[n.name] = meta_dic

    return meta


#####################################################
# date parsing and conversions
#####################################################
def parse_date(datein, fmt):
    from datetime import datetime
    import numpy as np
    try:
        if 'XX' in datein:
            min_date, max_date = ambiguous_date_to_date_range(datein, fmt)
            n_date = np.array((numerical_date(min_date), numerical_date(max_date)))
        else:
            tmp = datetime.strptime(datein, fmt).date()
            n_date = numerical_date(tmp)
    except:
        print("Can't parse ",datein)
        n_date=None

    return n_date


def numerical_date(date):
    from datetime import datetime
    days_in_year = date.toordinal()- datetime(year=date.year, month=1, day=1).date().toordinal()
    return date.year + days_in_year/365.25


def ambiguous_date_to_date_range(mydate, fmt):
    from datetime import datetime
    sep = fmt.split('%')[1][-1]
    min_date, max_date = {}, {}
    today = datetime.today().date()

    for val, field  in zip(mydate.split(sep), fmt.split(sep+'%')):
        f = 'year' if 'y' in field.lower() else ('day' if 'd' in field.lower() else 'month')
        if 'XX' in val:
            if f=='year':
                return None, None
            elif f=='month':
                min_date[f]=1
                max_date[f]=12
            elif f=='day':
                min_date[f]=1
                max_date[f]=31
        else:
            min_date[f]=int(val)
            max_date[f]=int(val)
    max_date['day'] = min(max_date['day'], 31 if max_date['month'] in [1,3,5,7,8,10,12]
                                           else 28 if max_date['month']==2 else 30)
    lower_bound = datetime(year=min_date['year'], month=min_date['month'], day=min_date['day']).date()
    upper_bound = datetime(year=max_date['year'], month=max_date['month'], day=max_date['day']).date()
    return (lower_bound, upper_bound if upper_bound<today else today)

##########################################
# IO
##########################################

def write_json(data, file_name, indent=1):
    import json
    try:
        handle = open(file_name, 'w')
    except IOError:
        pass
    else:
        json.dump(data, handle, indent=indent)
        handle.close()


########################################
# translation
#######################################3
nuc_alpha = 'ACGT-N'
aa_alpha = 'ACDEFGHIKLMNPQRSTVWY*-X'
TINY = 1e-12

def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'M'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Seq import Seq, CodonTable
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    from Bio.Data.CodonTable import TranslationError
    translation_exception = False

    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq(sequence).translate(gap='-'))
    except TranslationError:
        translation_exception = True

        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table  = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence)
        codons = np.fromstring(str_seq[:len(str_seq) - len(str_seq) % 3], dtype='S3')
        assert len(codons) > 0
        aas = []

        for c in codons:
            # Parse result of single codon translation, add amino acids as
            # appropriate.
            try:
                aa = codon_table.get(c)
                if aa is None:
                    if c == '---':
                        aas.append('-')
                    else:
                        aas.append('X')
                else:
                    aas.append(aa)
            except (TranslationError, ValueError):
                aas.append('X')

        translated_sequence = "".join(aas)

    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence


### get genes and translations
def get_genes_and_alignments(path, tree=True):
    import glob
    if tree:
        from filenames import tree_sequence_alignment as func
    else:
        from filenames import ref_alignment as func

    genes = []
    mask = func(path, prot='*')
    aln_files = glob.glob(mask)
    for aln_fname in aln_files:
        gene = aln_fname.rstrip(mask.split("*")[-1]).lstrip(mask.split('*')[0])
        genes.append((gene, aln_fname))
    return genes


def calc_af(aln_array, alpha):
    af = np.zeros((len(alpha), aln_array.shape[1]))
    for ai, state in enumerate(alpha):
        af[ai] += (aln_array==state).mean(axis=0)
    af[-1] = 1.0 - af[:-1].sum(axis=0)
    return af

def diversity_statistics(fname, nuc=True):
    from Bio import AlignIO
    aln_array = np.array(AlignIO.read(fname, 'fasta'))
    af = calc_af(aln_array, nuc_alpha if nuc else aa_alpha)
    tmp_af = af[:-2]/(af[:-2].sum(axis=0)+TINY)
    entropy = -(tmp_af*np.log(tmp_af+TINY)).sum(axis=0)

    return entropy


def load_features(reference, feature_names=None):
    from Bio import SeqIO
    features = {}
    for feat in SeqIO.read(reference, 'genbank').features:
        if feat.type=='CDS':
            if "locus_tag" in feat.qualifiers:
                fname = feat.qualifiers["locus_tag"][0]
                if feature_names is None or fname in feature_names:
                    features[fname] = feat

    return features

