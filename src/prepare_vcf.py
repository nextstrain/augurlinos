from __future__ import division, print_function
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import gzip
from shutil import copyfile
from filenames import (orig_meta_file_name, sequence_input, dropped_strains_file_name, recode_gzvcf_name,
    meta_file_name, ref_fasta, results_dir)
from util import write_sequence_meta_data, generic_argparse

def get_dropped_strains(path):
    fname = dropped_strains_file_name(path)
    dropped_strains = []
    if os.path.isfile(fname):
        with open(fname) as ifile:
            for line in ifile:
                fields = line.strip().split('#')
                if fields[0].strip():
                    dropped_strains.append(fields[0].strip())
    else:
        print("File with dropped strains not found. Looking for", fname)

    return dropped_strains


def mask_sites(path, ref, strip_loci):
    refSeq = SeqIO.parse(ref, format='fasta').next()
    refStr = str(refSeq.seq)

    maskArray = np.full(len(refStr),"0") #create masking Fasta for VCFTools to use
    maskedRef = np.array(list(refStr)) #to 'mask' the ref sequence

    df = pd.read_csv(strip_loci, sep='\t')
    for mi, m in df.iterrows():
        maskArray[(m.ChromStart-1):m.ChromEnd] = '1'
        maskedRef[(m.ChromStart-1):m.ChromEnd] = 'N'

    #need to get the CHROM name from the VCF file..
    #First try with normal open, then gzip - don't know file name at this point.
    #In new augur this can look at file ending - here cannot as fixed in filenames.py!!
    try:
        with gzip.open(recode_gzvcf_name(path)) as f:
            for line in f: #unknown quantity of comment at start of file
                if line[0] != '#':
                    line = line.strip()
                    header = line.split('\t')
                    chromName = header[0]
                    break
    except IOError as error:
        with open(recode_gzvcf_name(path)) as f:
            for line in f:
                if line[0] != '#':
                    line = line.strip()
                    header = line.split('\t')
                    chromName = header[0]
                    break

    exclude = []
    for i in xrange(len(maskArray)):
        if maskArray[i] == '1':
            exclude.append(chromName+"\t"+str(i+1))

    maskRefFile = ref+"_temp"
    with open(maskRefFile, 'w') as the_file:
        the_file.write("\n".join(exclude))

    maskedRef_seqRec = SeqRecord(seq=Seq("".join(maskedRef)),name=refSeq.name, id=refSeq.name, description='')

    #I stopped masking the Ref because I was unsure if the many columns of
    #N were causing problems later. Removing all variance and allowing these
    #regions to be the same as Ref should be the same, anyway.

    #with open(ref_fasta(path), "w") as output_handle:
    #    SeqIO.write(maskedRef_seqRec, output_handle, "fasta")

    return maskRefFile



if __name__ == '__main__':
    #to do - add so can pass vcf file instead of gzvcf file?

    import time
    start = time.time()

    parser = generic_argparse("parse vcf/vcf.gz file and meta_data to drop samples")
    parser.add_argument("--gzvcf", required=True, type=str,
                        help = "file with input sequences as gunzipped vcf")
    parser.add_argument("--ref", required=True, type=str,
                        help = "fasta file with reference sequence that vcf is mapped to")
    parser.add_argument("--strip_loci", required=False, type=str,
                        help = "file that contains loci to strip from analysis")
    args = parser.parse_args()
    path = args.path

    #set appropriate arg to VCFTools depending on input file
    #and set output params.
    #In new augur, will also want output ending (.vcf or .vcf.gz) set here
    vcfInCall = "--vcf"
    vcfOutCall = ""
    if args.gzvcf.endswith(('.gz', '.GZ')):
        vcfInCall = "--gzvcf"
        vcfOutCall = "| gzip -c"

    #VCFTools doesn't make directories if they don't exist
    if not os.path.isdir(results_dir(path)):
        os.makedirs(results_dir(path))

    #First copy to recode and rename and put in results
    #This is not *necessary* except to rename in augurlinos and move - but does not HAVE
    #to go through vcftools to be readable in further steps (build_tree.py) -- in new augur
    call = ["vcftools", vcfInCall, args.gzvcf, "--recode --stdout", vcfOutCall,">", recode_gzvcf_name(path)]
    print(" ".join(call))
    os.system(" ".join(call))

    ###See if want to mask sites
    if args.strip_loci:
        #This masks the ref sequence AND removes those sites from the VCF

        maskRefFile = mask_sites(path, args.ref, args.strip_loci)

        #for some reason vcftools doesn't seem to work if input and output are the same file
        #so create a copy we'll delete in a few lines
        #This uses the 'mask' file to remove the sites from the VCF
        copyfile(recode_gzvcf_name(path), recode_gzvcf_name(path)+"_temp")
        call = ["vcftools", "--exclude-positions", maskRefFile, vcfInCall, recode_gzvcf_name(path)+"_temp", "--recode --stdout", vcfOutCall, ">", recode_gzvcf_name(path)]
        print("Removing masked sites from VCF file... this may take some time.")
        os.system(" ".join(call))
        os.remove(recode_gzvcf_name(path)+"_temp")
        os.remove(maskRefFile)

    ###Get any sequences to be dropped
    dropped_strains = get_dropped_strains(path)

    #if some dropped, remove them
    if len(dropped_strains) != 0:
        #for some reason vcftools doesn't seem to work if input and output are the same file
        #so create a copy we'll delete in a few lines
        copyfile(recode_gzvcf_name(path), recode_gzvcf_name(path)+"_temp")
        toDrop = " ".join(["--remove-indv "+s for s in dropped_strains])
        call = ["vcftools", toDrop, vcfInCall, recode_gzvcf_name(path)+"_temp", "--recode --stdout", vcfOutCall, ">", recode_gzvcf_name(path)]
        print(" ".join(call))
        os.system(" ".join(call))
        os.remove(recode_gzvcf_name(path)+"_temp")

        #remove corresponding meta-data as well
        df = pd.read_csv(orig_meta_file_name(path), sep='\t')
        df = df[~df['strain'].isin(dropped_strains)]
        write_sequence_meta_data(path, df)


    #If the meta file hasn't already been written (b/c dropped strains), copy it
    if not os.path.isfile(meta_file_name(path)):
        copyfile(orig_meta_file_name(path), meta_file_name(path))

    #If the reference file hasn't already been written (b/c masking sites), copy it
    if not os.path.isfile(ref_fasta(path)):
        copyfile(args.ref, ref_fasta(path))

    os.remove('out.log') #remove vcftools log file thing.

    end = time.time()
    print("Prepare took {}".format(str(end-start)))

