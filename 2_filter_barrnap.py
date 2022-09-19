# Last updated march 10 2022

"""Filter the Barrnap outputs to keep the 16s sequences."""

import os
import sys
from Bio import SeqIO

"""Keeps just the best rRNAs sequences from a file."""
def filter_16s_seqs(fasta_file, taxon_name):

    bac_16s = []

    for seq_record in SeqIO.parse(fasta_file,'fasta'):
        # Checking that the sequences are 16S rRNAs.
        if seq_record.id[:2] == '16':

            seq_length = len(seq_record.seq)
            seq_cov = seq_record.id.split("_cov_")[1].split('.')[0]
            seq_number = seq_record.id.split('_')[2]
            seq_coords = seq_record.id.split(':')[-1]

            # Keep "high" expressed rRNAs that are longer.
            if int(seq_cov) >= 10:

                if int(seq_length) >= 900:

                    final_name = f'{taxon_name}_XX_Transcript_{seq_number}_' \
                        f'Len_{seq_length}_Cov_{seq_cov}_Pos_{seq_coords}'

                    bac_16s.append(f'>{final_name}\n{seq_record.seq}\n')
    # Save the filtered rRNAs for the next steps.
    with open(fasta_file.replace(".fas",".Filtered.fas"), "w+") as w:
        w.write("".join(bac_16s))

"""Runs the filtering step on many rRNA sequence files from Barrnap."""
def run_many_fasta(barrnap_folder):

    for fasta_file in os.listdir(barrnap_folder):

        if fasta_file.endswith('.Bac_16s_rRNA.fas'):

            taxon_name = fasta_file[:10]
            filter_16s_seqs(barrnap_folder+'/'+fasta_file, taxon_name)


if __name__ == '__main__':
    # Check to make sure that the barrnap folder is there.
    if len(sys.argv[1:]) == 1:
        barrnap_folder = sys.argv[1]
    else:
        print("\nUsage:\npython filter_barrnap.py [BARRNAP-FOLDER]\n")
        sys.exit()

    run_many_fasta(barrnap_folder)
