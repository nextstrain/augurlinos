def meta_file_name(path):
    return path+'/results/'+'meta.tsv'


def tree_meta_file_name(path):
    return path+'/results/'+'tree.tsv'


def sequence_input(path):
    return path+'/results/'+'orig.fasta'


def raw_alignment(path):
    return path+'/results/'+'nuc_aln.fasta'


def ref_alignment(path, prot='nuc'):
    return path+'/results/'+'ref_%s_aln.fasta'%prot


def tree_sequence_alignment(path, prot='nuc'):
    return path+'/results/'+'tree_%s_aln.fasta'%prot


def tree_newick(path):
    return path+'/results/'+'tree.nwk'


def sequence_gtr_model(path):
    return path+'/results/'+'seq_gtr.txt'


def mugration_model(path, field):
    return path+'/results/'+'%s_gtr.txt'%field



#######################################
# Auspice json file name
#######################################

def sequence_json(path,prefix):
    return path+'/auspice/'+prefix+'_sequence.json'

def tree_json(path, prefix):
    return path+'/auspice/'+prefix+'_tree.json'

def diversity_json(path, prefix):
    return path+'/auspice/'+prefix+'_entropy.json'
