######################################
# input files
######################################

def dropped_strains_file_name(path):
    return path+'/data/'+'dropped_strains.txt'


def color_maps(path):
    return path+'/data/'+'color.tsv'


#EBH 29 Nov 17
def orig_meta_file_name(path):
    return path+'/data/'+'meta.tsv'

######################################
# process files
######################################

#EBH 29 Nov 17
def recode_gzvcf_name(path):
    return path+'/results/'+'seqs.vcf.gz'


#EBH 29 Nov 17
def ref_fasta(path):
    return path+'/results/'+'ref.fasta'
    

#EBH 30 Nov 17
def var_site_alignment(path):
    return path+'/results/'+'var_sites.fasta'
    
    
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

    
#EBH 7 Dec 17
def tree_vcf_alignment(path, prot='nuc'):
    return path+'/results/'+'tree_%s_aln.vcf'%prot
    

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
    return path+'/auspice/'+prefix+'_sequences.json'

def tree_json(path, prefix):
    return path+'/auspice/'+prefix+'_tree.json'

def diversity_json(path, prefix):
    return path+'/auspice/'+prefix+'_entropy.json'
