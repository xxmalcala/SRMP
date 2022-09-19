
# Last updated march 02 2022

"""Runs barrnap to identify 16s rRNA sequences from a FOLDER of FASTA files."""

import os
import sys

# Runs the barrnap program to find the bacterial rRNAs.
def run_barrnap(fasta_file, unique_out_name):
    barrnap_cmd = 'barrnap --kingdom bac ' \
        '--threads 28 ' \
        '--lencutoff 0.6 ' \
        '--outseq '+unique_out_name+'.Bac_16s_rRNA.fas ' +\
        fasta_file
    # Runs the barrnap progam with the command.
    os.system(barrnap_cmd)

# Runs Barrnap on all the transcriptomes in a folder.
def run_many_fasta(fasta_folder):

    for fasta_file in os.listdir(fasta_folder):

        if fasta_file.endswith('.fasta'):
            # We use abbreviations for the organism to make it easier to track.
            taxon_name = fasta_file[:10]
            run_barrnap(fasta_folder+'/'+fasta_file, taxon_name)

if __name__ == '__main__':
    if len(sys.argv[1:]) == 1:
        fasta_folder = sys.argv[1]
    else:
        print("\nUsage:\npython example_barrnap.py [FASTA-FOLDER]\n")
        sys.exit()

    run_many_fasta(fasta_folder)
