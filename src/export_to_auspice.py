import numpy as np
from util import (read_sequence_meta_data, read_tree_meta_data,
                  get_genes_and_alignments, generic_argparse)
from filenames import tree_newick, tree_json, tree_sequence_alignment, sequence_json, diversity_json
from util import write_json, load_features, diversity_statistics
from Bio import Phylo


def tree_to_json(node, extra_attr = []):
    tree_json = {}
    str_attr = ['strain']
    num_attr = ['num_date']
    if hasattr(node, 'name'):
        tree_json['strain'] = node.name

    for prop in str_attr:
        if hasattr(node, prop):
            tree_json[prop] = node.__getattribute__(prop)
    for prop in num_attr:
        if hasattr(node, prop):
            try:
                tree_json[prop] = round(node.__getattribute__(prop),5)
            except:
                print("cannot round:", node.__getattribute__(prop), "assigned as is")
                tree_json[prop] = node.__getattribute__(prop)

    for prop in extra_attr:
        if len(prop)==2 and callable(prop[1]):
            if hasattr(node, prop[0]):
                tree_json[prop] = prop[1](node.__getattribute__(prop[0]))
        else:
            if hasattr(node, prop):
                tree_json[prop] = node.__getattribute__(prop)

    tree_json['tvalue'] = tree_json['num_date']
    if node.clades:
        tree_json["children"] = []
        for ch in node.clades:
            tree_json["children"].append(tree_to_json(ch, extra_attr))
    return tree_json


def attach_tree_meta_data(T, node_meta):
    def parse_mutations(muts):
        return muts.split(',') if type(muts)==str else ""

    for n in T.find_clades(order='preorder'):
        n.attr={}
        n.aa_muts={}
        for field, val in node_meta[n.name].items():
            if 'mutations' in field:
                if field=='mutations':
                    muts = parse_mutations(val)
                    if muts:
                        n.__setattr__('muts', muts)
                else:
                    prot = '_'.join(field.split('_')[:-1])
                    muts = parse_mutations(val)
                    if muts:
                        n.aa_muts[prot] = muts
            elif field in ['branch_length', 'mutation_length', 'clock_length',
                           'clade', 'num_date']:
                n.__setattr__(field, val)
                n.attr[field] = val
            else:
                n.attr[field] = val

    T.root.attr['div']=0
    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            bl =  n.mutation_length if hasattr(n, "mutation_length") else "branch_length"
            c.attr["div"] = n.attr["div"] + bl



def export_sequence_json(T, path, prefix, indent=None):
    from Bio import SeqIO
    plain_export = 0.99

    elems = {'root':{}}
    for node in T.find_clades():
        elems[node.clade] = {}

    for gene, aln_fname in get_genes_and_alignments(path, tree=True):
        seqs={}
        for seq in SeqIO.parse(aln_fname, 'fasta'):
            seqs[seq.name] = seq

        root_seq = seqs[T.root.name]
        elems['root'][gene] = "".join(root_seq)
        for node in T.find_clades():
            nseq = seqs[node.name]
            if hasattr(node, "clade"):
                differences = {pos:state for pos, (state, ancstate) in
                            enumerate(zip(nseq, elems['root'][gene]))
                            if state!=ancstate}
                if len(differences)<=plain_export*len(seq):
                    elems[node.clade][gene] = differences
                else:
                    elems[node.clade][gene] = seq

    fname = sequence_json(path, prefix)
    write_json(elems, fname, indent=indent)


def export_metadata_json(T, path, prefix, indent):
    print("Writing out metaprocess")
    meta_json = {}

    meta_json["virus_count"] = T.count_terminals()
    from datetime.date import today
    meta_json["updated"] = today().strftime('%Y-%m-%d')
    meta_json["author_info"] = {}
    meta_json["seq_author_map"] = {}


    # join up config color options with those in the input JSONs.
    col_opts = process.config["auspice"]["color_options"]
    if process.colors:
        for trait, col in process.colors.iteritems():
            if trait in col_opts:
                col_opts[trait]["color_map"] = col
            else:
                process.log.warn("{} in colors (input JSON) but not auspice/color_options. Ignoring".format(trait))

    meta_json["color_options"] = col_opts
    if "date_range" in process.config["auspice"]:
        meta_json["date_range"] = process.config["auspice"]["date_range"]
    if "analysisSlider" in process.config["auspice"]:
        meta_json["analysisSlider"] = process.config["auspice"]["analysisSlider"]
    meta_json["panels"] = process.config["auspice"]["panels"]
    meta_json["updated"] = time.strftime("X%d %b %Y").replace('X0','X').replace('X','')
    meta_json["title"] = process.info["title"]
    meta_json["maintainer"] = process.info["maintainer"]
    meta_json["filters"] = process.info["auspice_filters"]

    if "defaults" in process.config["auspice"]:
        meta_json["defaults"] = process.config["auspice"]["defaults"]

    try:
        from pygit2 import Repository, discover_repository
        current_working_directory = os.getcwd()
        repository_path = discover_repository(current_working_directory)
        repo = Repository(repository_path)
        commit_id = repo[repo.head.target].id
        meta_json["commit"] = str(commit_id)
    except ImportError:
        meta_json["commit"] = "unknown"
    if len(process.config["auspice"]["controls"]):
        meta_json["controls"] = process.make_control_json(process.config["auspice"]["controls"])
    meta_json["geo"] = process.lat_longs
    write_json(meta_json, prefix+'_meta.json')


def export_diversity(path, prefix, reference, indent=None):
    '''
    write the alignment entropy of each alignment (nucleotide and translations) to file
    '''
    genes = load_features(reference)
    entropy_json = {}
    for feat, aln_fname in get_genes_and_alignments(path, tree=False):
        entropy = diversity_statistics(aln_fname, nuc=feat=='nuc')
        S = [max(0,round(x,4)) for x in entropy]
        n = len(S)
        if feat=='nuc':
            entropy_json[feat] = {'pos':range(0,n), 'codon':[x//3 for x in range(0,n)], 'val':S}
        elif feat in genes:
            entropy_json[feat] = {'pos':[x for x in genes[feat]][::3],
                                  'codon':range(n), 'val':S}
    write_json(entropy_json, diversity_json(path, prefix), indent=indent)


def tree_layout(T):
    yval=T.count_terminals()
    for n in T.find_clades(order='postorder'):
        if n.is_terminal():
            n.yvalue=yval
            yval-=1
        else:
            child_yvalues = [c.yvalue for c in n]
            n.yvalue=0.5*(np.min(child_yvalues)+np.max(child_yvalues))
        n.xvalue = n.attr['div']



if __name__ == '__main__':
    parser =  generic_argparse("Export precomputed data as auspice jsons")
    parser.add_argument('--prefix', required=True,
                        help="prefix for json files that are passed on to auspice (e.g., zika.fasta)")
    parser.add_argument('--reference', required=True,
                        help="reference sequence needed for entropy feature export")

    args = parser.parse_args()
    path = args.path

    T = Phylo.read(tree_newick(path), 'newick')
    seq_meta = read_sequence_meta_data(path)
    tree_meta = read_tree_meta_data(path)
    attach_tree_meta_data(T, tree_meta)
    tree_layout(T)
    fields_to_export = tree_meta.values()[0].keys()+["tvalue","yvalue", "xvalue", "attr","muts", "aa_muts"]
    tjson = tree_to_json(T.root, extra_attr=fields_to_export)
    write_json(tjson, tree_json(path, args.prefix))

    export_sequence_json(T, path, args.prefix, indent=1)

    export_diversity(path, args.prefix, args.reference, indent=1)
