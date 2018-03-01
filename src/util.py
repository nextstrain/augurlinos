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
            seq.seq.upper()
            SeqIO.write(seq, ofile, format='fasta')


def read_sequence_meta_data(path):
    from filenames import meta_file_name
    import pandas as pd
    df = pd.read_csv(meta_file_name(path), sep='\t').fillna('')
    #import ipdb; ipdb.set_trace()
    return {m[0]:m.to_dict() for mi, m in df.iterrows()}


def write_sequence_meta_data(path, df):
    from filenames import meta_file_name
    df.to_csv(meta_file_name(path), sep='\t', index=False)


def read_tree_meta_data(path):
    from filenames import tree_meta_file_name
    # import pandas as pd
    # df = pd.read_csv(tree_meta_file_name(path), sep='\t')
    # return {m[0]:m.to_dict() for mi, m in df.iterrows()}
    import json
    with open(tree_meta_file_name(path)) as ifile:
        meta = json.load(ifile)
    return meta


def write_tree_meta_data(path, meta_dic, indent=1):
    from filenames import tree_meta_file_name
    # import pandas as pd
    # header = set()
    # for m in meta_dic.values():
    #     header.update(m.keys())
    # header = sorted(header, key=lambda x:x!='name')
    # df = pd.DataFrame(meta_dic.values(), columns = header)
    # df.to_csv(tree_meta_file_name(path), sep='\t', index=False)
    import json
    with open(tree_meta_file_name(path), 'w') as ofile:
        json.dump(meta_dic, ofile, indent=indent)



def collect_tree_meta_data(T, fields, isvcf=False, meta=None):
    def mutation_format(muts):
        if isvcf:
            #converts from Python numbering [0] to standard [1] for output
            return ",".join(['%s%d%s'%(x[0], x[1]+1, x[2]) for x in muts])
        else:
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

def read_in_vcf(vcf_file, ref_file, compressed=True):
    """
    Reads in a vcf.gz file (or vcf if compressed is False) and associated
    reference sequence fasta (to which the VCF file is mapped)

    Parses mutations, insertions, and deletions and stores them in a nested dict
    with the format:
    {'reference':'AGCTCGA..A',
     'sequences': { 'seq1':{4:'A', 7:'-'}, 'seq2':{100:'C'} },
     'insertions': { 'seq1':{4:'ATT'}, 'seq3':{1:'TT', 10:'CAG'} },
     'positions': [1,4,7,10,100...] }

    Calls with values 0/1 (or 0/2, etc) are ignored.
    Positions are stored to correspond the location in the reference sequence
    in Python (numbering is transformed to start at 0)

    Args
    ----
    vcf_file : string
        Path to the vcf or vcf.gz file to be read in
    ref_file : string
        Path to the fasta reference file to be read in
    compressed : boolean
        Specify false if VCF file is not compressed (not vcf.gz)

    Returns
    --------
    compress_seq : nested dict
        Contains the following keys:

        references : string
            String of the reference sequence read from the Fasta
        sequences : nested dict
            Dict containing sequence names as keys which map to dicts
            that have position as key and the single-base mutation (or deletion)
            as values
        insertions : nested dict
            Dict in the same format as the above, which stores insertions and their
            locations. The first base of the insertion is the same as whatever is
            currently in that position (Ref if no mutation, mutation in 'sequences'
            otherwise), so the current base can be replaced by the bases held here
            without losing that base.
        positions : list
            Python list of all positions with a mutation, insertion, or deletion.

    EBH 4 Dec 2017
    """
    #vcf is inefficient for the data we want. House code is *much* faster.
    import gzip
    from Bio import SeqIO
    import numpy as np

    sequences = {}
    insertions = {}
    positions = []

    nsamp = 0
    posLoc = 0
    refLoc = 0
    altLoc = 0
    sampLoc = 9

    if compressed: #must use 2 diff functions depending on compressed or not
        opn = gzip.open
    else:
        opn = open

    #with gzip.open(vcf_file) as f:
    with opn(vcf_file) as f:
        for line in f:
            if line[0] != '#':
                #actual data - most common so first in 'if-list'!
                line = line.strip()
                dat = line.split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                calls = np.array(dat[sampLoc:])


                #get samples that differ from Ref at this site
                recCalls = {}
                k=0
                for sa in calls:
                    if ':' in sa: #if proper VCF file
                        gt = sa.split(':')[0]
                    else: #if 'pseudo' VCF file (nextstrain output)
                        gt = sa
                    if gt != '.' and gt[0] != '.' and gt[0] != '0':
                        recCalls[samps[k]] = gt
                    k+=1

                #store the position and the alt
                for seq, gen in recCalls.iteritems():
                    if gen[0] != '0' and gen[2] != '0':
                        #if is 0/1 or 1/0, ignore - uncertain call
                        alt = str(ALT[int(gen[0])-1])   #get the index of the alternate
                        ref = REF
                        pos = POS-1     #VCF numbering starts from 1, but Reference seq numbering
                                        #will be from 0 because it's python!

                        if seq not in sequences.keys():
                            sequences[seq] = {}

                        #figure out if insertion or deletion
                        #insertion where there is also deletions (special handling)
                        if len(ref) > 1 and len(alt)>len(ref):
                            #print "nonstandard insertion at pos {}".format(record.POS)
                            if seq not in insertions.keys():
                                insertions[seq] = {}
                            for i in xrange(len(ref)):
                                #if the pos doesn't match, store in sequences
                                if ref[i] != alt[i]:
                                    sequences[seq][pos+i] = alt[i]
                                    #if pos+1 not in positions:
                                    positions.append(pos+i)
                                #if about to run out of ref, store rest:
                                if (i+1) >= len(ref):
                                    insertions[seq][pos+i] = alt[i:]
                                    #print "at pos {}, storing {} at pos {}".format(record.POS, alt[i:], (pos+i))

                        #deletion
                        elif len(ref) > 1:
                            for i in xrange(len(ref)):
                                #if ref is longer than alt, these are deletion positions
                                if i+1 > len(alt):
                                    sequences[seq][pos+i] = '-'
                                    #if pos+i not in positions:
                                    positions.append(pos+i)
                                #if not, there may be mutations
                                else:
                                    if ref[i] != alt[i]:
                                        if alt[i] == '.':
                                            sequences[seq][pos+i] = '-'
                                        else:
                                            sequences[seq][pos+i] = alt[i]
                                        #if pos+i not in positions:
                                        positions.append(pos+i)

                        #insertion
                        elif len(alt) > 1:
                            #keep a record of insertions so can put them in if we want, later
                            if seq not in insertions.keys():
                                insertions[seq] = {}
                            insertions[seq][pos] = alt
                            #First base of insertions always matches ref, so don't need to store

                        #no indel
                        else:
                            sequences[seq][pos] = alt
                            #if pos not in positions:
                            positions.append(pos)


            elif line[0] == '#' and line[1] == 'C':
                #header line, get all the information
                line = line.strip()
                header = line.split('\t')
                headNP = np.array(header)
                posLoc = np.where(headNP=='POS')[0][0]
                refLoc = np.where(headNP=='REF')[0][0]
                altLoc = np.where(headNP=='ALT')[0][0]
                sampLoc = np.where(headNP=='FORMAT')[0][0]+1
                samps = header[sampLoc:]
                nsamp = len(samps)


            #else you are a comment line, ignore.



    positions = np.array(positions)
    positions = np.unique(positions)
    positions = np.sort(positions)

    refSeq = SeqIO.parse(ref_file, format='fasta').next()
    refSeqStr = str(refSeq.seq)

    compress_seq = {'reference':refSeqStr,
                    'sequences': sequences,
                    'insertions': insertions,
                    'positions': positions}

    return compress_seq

def read_in_translate_vcf(vcf_file, ref_file, compressed=True):
    """
    Reads in a vcf file where TRANSLATIONS have been stored and associated
    reference sequence fasta (to which the VCF file is mapped)
    This is the file output by "write_VCF_translation" below

    Returns a nested dict in the same format as is *input* in "write_VCF_translation" below,
    with a nested dict for each gene, which contains 'sequences', 'positions', and 'reference'
    """
    from Bio import SeqIO
    import numpy as np

    prots = {}

    posLoc = 0
    refLoc = 0
    altLoc = 0
    sampLoc = 9

    with open(vcf_file) as f:
        for line in f:
            if line[0] != '#':
                #actual data
                line = line.strip()
                dat = line.split('\t')
                POS = int(dat[posLoc])
                REF = dat[refLoc]
                ALT = dat[altLoc].split(',')
                GEN = dat[0] #'CHROM' or the gene name here
                calls = np.array(dat[sampLoc:])

                #get samples that differ from Ref at this site
                recCalls = {}
                k=0
                for sa in calls:
                    if sa != '.':
                        recCalls[samps[k]] = sa
                    k+=1

                #store position and the altLoc
                for seq, gen in recCalls.iteritems():
                    alt = str(ALT[int(gen[0])-1])   #get the index of the alternate
                    ref = REF
                    pos = POS-1     #VCF numbering starts from 1, but Reference seq numbering
                                    #will be from 0 because it's python!
                    gen = GEN       #from CHROM, gene name

                    if gen not in prots.keys():
                        prots[gen] = {}
                        prots[gen]['sequences'] = {}
                        prots[gen]['positions'] = []
                        prots[gen]['reference'] = ''
                    if seq not in prots[gen]['sequences'].keys():
                        prots[gen]['sequences'][seq] = {}

                    #will never be insertion or deletion! because translation.
                    prots[gen]['sequences'][seq][pos] = alt
                    prots[gen]['positions'].append(pos)

            elif line[0] == '#' and line[1] == 'C':
                #header line, get all the information
                line = line.strip()
                header = line.split('\t')
                headNP = np.array(header)
                posLoc = np.where(headNP=='POS')[0][0]
                refLoc = np.where(headNP=='REF')[0][0]
                altLoc = np.where(headNP=='ALT')[0][0]
                sampLoc = np.where(headNP=='FORMAT')[0][0]+1 #first sample should be 1 after FORMAT column
                samps = header[sampLoc:]
                nsamp = len(samps)


    for refSeq in SeqIO.parse(translation_ref_file(path), format='fasta'):
        prots[refSeq.name]['reference'] = str(refSeq.seq)
        posN = np.array(prots[refSeq.name]['positions'])
        posN = np.unique(posN)
        prots[refSeq.name]['positions'] = np.sort(posN)

    return prots


def read_in_DRMs(drm_file):
    """
    Reads in a tab-delim file giving information on Drug Resistance Mutations.
    Format at the moment must include columns titled:
    GENOMIC_POSITION (bp of the mutation)
    ALT_BASE    (the alt-base conferring resistance)
    SUBSTITUTION    (the aminoacid-codonloc-aminoacid change)
    DRUG    (the drug of resistance)
    GENE    (the gene the mutation is in)

    Other columns can be present but will be ignored

    Returns a nested dict containing:
        drmPositions: an array of all the positions of DRMs

        DRMs: another nested dict with information about the mutations
            'position' is the key, and maps to a nested dict that contains:
            base, drug, AA, gene
            base and AA are arrays as one bp position can have multiple 'alt' bases

    """
    import pandas as pd
    import numpy as np

    DRMs = {}
    drmPositions = []

    df = pd.read_csv(drm_file, sep='\t')
    for mi, m in df.iterrows():
        pos = m.GENOMIC_POSITION-1 #put in python numbering
        drmPositions.append(pos)

        if pos in DRMs:
            DRMs[pos]['base'].append(m.ALT_BASE)
            DRMs[pos]['AA'].append(m.SUBSTITUTION)
        else:
            DRMs[pos] = {}
            DRMs[pos]['base'] = [m.ALT_BASE]
            DRMs[pos]['drug'] = m.DRUG
            DRMs[pos]['AA'] = [m.SUBSTITUTION]
            DRMs[pos]['gene'] = m.GENE

    drmPositions = np.array(drmPositions)
    drmPositions = np.unique(drmPositions)
    drmPositions = np.sort(drmPositions)

    DRM_info = {'DRMs': DRMs,
            'drmPositions': drmPositions}

    return DRM_info


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
    import os

    #in case auspice folder does not exist yet
    if not os.path.exists(os.path.dirname(file_name)):
        try:
            os.makedirs(os.path.dirname(file_name))
        except OSError: #Guard against race condition
            if not os.path.isdir(os.path.dirname(file_name)):
                raise
    try:
        handle = open(file_name, 'w')
    except IOError:
        raise
    else:
        json.dump(data, handle, indent=indent)
        handle.close()

def write_VCF_style_alignment(tree_dict, file_name):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the alignment from the input of a dict
    in a similar format to what's created from the read_in_vcf function above.

    EBH 7 Dec 2017
    """
    sequences = tree_dict['sequences']
    ref = tree_dict['reference']
    positions = tree_dict['positions']

    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+sequences.keys()
    with open(file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    vcfWrite = []

    #now get the variable positions and calls for every sample (node)
    i=0
    while i < len(positions):
        #get the pattern at this position
        #but only if not deletion!
        #VCF handles deletions in a really weird way, so we have to accomodate this...
        pi = positions[i]
        pos = pi+1 #change numbering to match VCF not python
        refb = ref[pi] #reference base at this position
        delete = False #deletion at this pos, need to grab previous pos too (which is invariable)
        deleteGroup = False #deletion at next pos (mutation at this pos)

        #I've updated how getting patterns is handled - it turns out using an 'if'
        #is not v. efficient, better to try every location and handle exceptions,
        #as on average a variable location will not be variable for one sequence.

        pattern = []
        pattern2 = [] #go ahead and get both, so only loop once!
        for k,v in sequences.iteritems():
            try:
                pattern.append(sequences[k][pi])
            except KeyError, e:
                pattern.append(ref[pi])

            try:
                pattern2.append(sequences[k][pi+1])
            except KeyError, e:
                pattern2.append(ref[pi+1])

        pattern = np.array(pattern)
        pattern2 = np.array(pattern2)

        #pattern = np.array([ sequences[k][pi] if pi in sequences[k].keys() else '.' for k,v in sequences.iteritems() ])  #old way

        #if a deletion here, need to gather up all bases, and position before
        if any(pattern == '-'):
            if pos != 1:
                deleteGroup = True
                delete = True
            else:
                #if theres a deletion in 1st pos, VCF files do not handle this well.
                #proceed keeping it as '-' for alt, but warn user to check output.
                print "WARNING: You have a deletion in the first position of your alignment. VCF format does not handle this well. Please check the output to ensure it is correct."
        else:
            #if there's a deletion in next pos, need gather up bases
            #pattern2 = np.array([ sequences[k][pi+1] if pi+1 in sequences[k].keys() else '.' for k,v in sequences.iteritems() ])
            if any(pattern2 == '-'):
                deleteGroup = True

        #if there is a deletion, treat affected bases as 1 'call':
        if delete or deleteGroup:
            if delete: #need to get the position before!
                i-=1
                pi-=1
                pos = pi+1
                refb = ref[pi]
                #re-get pattern if this is the case, same new method as above:
                pattern = []
                for k,v in sequences.iteritems():
                    try:
                        pattern.append(sequences[k][pi])
                    except KeyError, e:
                        pattern.append(ref[pi])
                pattern = np.array(pattern)

            sites = []
            sites.append(pattern)

            #gather all positions affected by deletion
            while positions[i+1] == pi+1:
                i+=1
                pi = positions[i]
                #again, new method
                pattern = []
                for k,v in sequences.iteritems():
                    try:
                        pattern.append(sequences[k][pi])
                    except KeyError, e:
                        pattern.append(ref[pi])
                pattern = np.array(pattern)
                #pattern = np.array([ sequences[k][pi] if pi in sequences[k].keys() else ref[pi] for k,v in sequences.iteritems() ])
                if any(pattern == '-'): #if part of deletion, append
                    sites.append(pattern)
                    refb = refb+ref[pi]
                else: #this is another mutation next to the deletion!
                    i-=1    #don't append, leave for next loop!

            #group them into 'calls'
            sites = np.asarray(sites)
            align = np.rot90(sites)
            align = np.flipud(align)

            #get rid of deletions, and put '.' for calls that match ref
            fullpat = []
            for pt in align:
                pat = "".join(pt).replace('-','')
                if pat == refb:
                    fullpat.append('.')
                else:
                    fullpat.append(pat)

            pattern = np.array(fullpat)

        else: #if no deletion, need to replace ref with '.' as in VCF files
            pattern[pattern==refb] = '.'

        #get the list of ALTs - minus any '.'!
        uniques = np.unique(pattern)
        uniques = uniques[np.where(uniques!='.')]

        #Convert bases to the number that matches the ALT
        j=1
        for u in uniques:
            pattern[np.where(pattern==u)[0]] = str(j)
            j+=1
        #Now convert these calls to #/# (VCF format)
        calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]
        if len(uniques)==0:
            print "util.py: UNEXPECTED ERROR WHILE CONVERTING TO VCF AT POSITION {}. No variation!".format(str(pi))
            #break   #keep going, as every other site should be fine...

        #put it all together and write it out!
        #increment positions by 1 so it's in VCF numbering not python numbering
        output = ["MTB_anc", str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls

        vcfWrite.append("\t".join(output))

        #with open(file_name, 'a') as the_file:
        #    the_file.write("\t".join(output)+"\n")
        i+=1

    with open(file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))


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
    import os.path
    from filenames import translation_ref_file, translation_vcf_file

    genes = []

    if os.path.isfile(translation_ref_file(path)): #vcf file, use different method!
        from Bio import SeqIO
        for seq in SeqIO.parse(translation_ref_file(path), format='fasta'):
            genes.append((seq.name, translation_vcf_file(path)))

    else:
        import glob
        if tree:
            from filenames import tree_sequence_alignment as func
        else:
            from filenames import ref_alignment as func

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
    #read in appropriately whether GFF or Genbank
    #checks explicitly for GFF otherwise assumes Genbank
    features = {}

    if '.gff' in reference.lower():
        #looks for 'gene' and 'gene' as best for TB
        from BCBio import GFF
        limit_info = dict( gff_type = ['gene'] )

        in_handle = open(reference)
        for rec in GFF.parse(in_handle, limit_info=limit_info):
            for feat in rec.features:
                if "gene" in feat.qualifiers:
                    fname = feat.qualifiers["gene"][0]
                else:
                    fname = feat.qualifiers["locus_tag"][0]
                if feature_names is None or fname in feature_names:
                    features[fname] = feat

        if feature_names is not None:
            for fe in feature_names:
                if fe not in features:
                    print "Couldn't find gene {} in GFF or GenBank file".format(fe)

        in_handle.close()

    else:
        from Bio import SeqIO
        for feat in SeqIO.read(reference, 'genbank').features:
            if feat.type=='CDS':
                if "locus_tag" in feat.qualifiers:
                    fname = feat.qualifiers["locus_tag"][0]
                    if feature_names is None or fname in feature_names:
                        features[fname] = feat

    return features


def write_VCF_translation(prot_dict, vcf_file_name, ref_file_name):
    """
    Writes out a VCF-style file (which seems to be minimally handleable
    by vcftools and pyvcf) of the AA differences between sequences and the reference.
    This is a similar format created/used by read_in_vcf except that there is one
    of these dicts (with sequences, reference, positions) for EACH gene.

    Also writes out a fasta of the reference alignment.

    EBH 12 Dec 2017
    """

    #for the header
    seqNames = prot_dict[prot_dict.keys()[0]]['sequences'].keys()

    #prepare the header of the VCF & write out
    header=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+seqNames
    with open(vcf_file_name, 'w') as the_file:
        the_file.write( "##fileformat=VCFv4.2\n"+
                        "##source=NextStrain_Protein_Translation\n"+
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        the_file.write("\t".join(header)+"\n")

    refWrite = []
    vcfWrite = []

    #go through for every gene/protein
    for fname, prot in prot_dict.iteritems():
        sequences = prot['sequences']
        ref = prot['reference']
        positions = prot['positions']

        #write out the reference fasta
        refWrite.append(">"+fname)
        refWrite.append(ref)

        #go through every variable position
        #There are no deletions here, so it's simpler than for VCF nuc sequenes!
        i=0
        while i < len(positions):
            pi = positions[i]
            pos = pi+1 #change numbering to match VCF not python
            refb = ref[pi] #reference base at this position

            pattern = np.array([ sequences[k][pi] if pi in sequences[k].keys() else '.' for k,v in sequences.iteritems() ])

            #get the list of ALTs - minus any '.'!
            uniques = np.unique(pattern)
            uniques = uniques[np.where(uniques!='.')]

            #Convert bases to the number that matches the ALT
            j=1
            for u in uniques:
                pattern[np.where(pattern==u)[0]] = str(j)
                j+=1
            #Now convert these calls to #/# (VCF format)
            calls = [ j+"/"+j if j!='.' else '.' for j in pattern ]
            if len(uniques)==0:
                print "UNEXPECTED ERROR WHILE CONVERTING TO VCF AT POSITION {}".format(str(pi))
                break

            #put it all together and write it out!
            #increment positions by 1 so it's in VCF numbering not python numbering
            output = [fname, str(pos), ".", refb, ",".join(uniques), ".", "PASS", ".", "GT"] + calls

            vcfWrite.append("\t".join(output))

            i+=1

    #write it all out
    with open(ref_file_name, 'w') as the_file:
        the_file.write("\n".join(refWrite))

    with open(vcf_file_name, 'a') as the_file:
        the_file.write("\n".join(vcfWrite))


def load_lat_long_defs():
    places = {}
    with open('../fauna/source-data/geo_lat_long.tsv', 'r') as latlongfile:
        header = latlongfile.readline().strip().split('\t')
        for line in latlongfile:
            try:
                place, country, latitude, longitude = line.strip().split('\t')
                places[place] = {'country_code':country,
                                 'latitude':float(latitude),
                                 'longitude':float(longitude)}
            except:
                print("trouble parsing", line)
    return places