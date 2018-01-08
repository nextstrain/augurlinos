from __future__ import division, print_function
import os
from Bio import SeqIO
import pandas as pd
from filenames import orig_meta_file_name, sequence_input, dropped_strains_file_name, recode_gzvcf_name, meta_file_name, ref_fasta
from util import write_sequence_meta_data, generic_argparse
from shutil import copyfile

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



if __name__ == '__main__':
    #to do - add so can pass vcf file instead of gzvcf file?

    import time
    start = time.time()

    parser = generic_argparse("parse gzvcf file and meta_data to drop samples")
    parser.add_argument("--gzvcf", required=True, type=str,
                        help = "file with input sequences as gunzipped vcf")
    parser.add_argument("--ref", required=True, type=str,
                        help = "fasta file with reference sequence that vcf is mapped to")
    args = parser.parse_args()
    path = args.path

    #First copy to recode and rename and put in results
    call = ["vcftools --gzvcf", args.gzvcf, "--recode --stdout | gzip -c >", recode_gzvcf_name(path)]
    print("VCFTools call:")
    print(" ".join(call))
    os.system(" ".join(call))

    #Get any sequences to be dropped
    dropped_strains = get_dropped_strains(path)

    #if some dropped, remove them in a loop
    if len(dropped_strains) != 0:
        #for some reason vcftools doesn't seem to work if input and output are the same file
        #so create a copy we'll delete in a few lines
        copyfile(recode_gzvcf_name(path), recode_gzvcf_name(path)+"_temp")
        toDrop = " ".join(["--remove-indv "+s for s in dropped_strains])
        call = ["vcftools", toDrop, "--gzvcf", recode_gzvcf_name(path)+"_temp", "--recode --stdout | gzip -c >", recode_gzvcf_name(path)]
        print(" ".join(call))
        os.system(" ".join(call))
        os.remove(recode_gzvcf_name(path)+"_temp")

        #remove corresponding meta-data as well
        df = pd.read_csv(orig_meta_file_name(path), sep='\t')
        df = df[~df['strain'].isin(dropped_strains)]
        write_sequence_meta_data(path, df)

    else:
        #if none dropped, just copy the meta file to results
        copyfile(orig_meta_file_name(path), meta_file_name(path))


    #finally copy the reference fasta so we know where it is and what its called
    copyfile(args.ref, ref_fasta(path))

    end = time.time()
    print("Prepare took {}".format(str(end-start)))

