# Last updated march 16 2022

"""Merge the sequences from the filtered Barrnap data. This keeps just the
unique rRNAs."""

import os
import sys


def merge_fasta_files(filtered_folder):
    # Merges all the filtered rRNA files together into a single file called
    # Merged.Bac_16s_rRNA.Filtered.fas

    os.system(f'cat {filtered_folder}/*Filtered.fas > ' \
        f'{filtered_folder}/Merged.Bac_16s_rRNA.Filtered.fas')

    return f'{filtered_folder}/Merged.Bac_16s_rRNA.Filtered.fas'

"""Step keeps just the unique rRNA sequences. We allow a 1% error, just in case
of sequencing mistakes."""
def cluster_merged_fasta(merged_fasta, id = 0.99, threads = 4):

    perc_id = int(id*100)

    output_file_name = f'Merged.16s_rRNA.{perc_id}perc_identity.fas'

    cd_hit_cmd = f'cd-hit-est -G 0 ' \
                f'-c {id} ' \
                f'-d 0 ' \
                f'-aS 1.0 ' \
                f'-aL 0.005 ' \
                f'-T {threads} ' \
                f'-i {merged_fasta} ' \
                f'-o {output_file_name}'

    os.system(cd_hit_cmd)


if __name__ == '__main__':
    if len(sys.argv[1:]) == 1:
        filtered_folder = sys.argv[1]
    else:
        print("\nUsage:\npython cluster_16s_rRNA.py [BARRNAP-FOLDER]\n")
        sys.exit()

    merged_fasta = merge_fasta_files(filtered_folder)
    cluster_merged_fasta(merged_fasta, id = 0.99, threads = 4)
