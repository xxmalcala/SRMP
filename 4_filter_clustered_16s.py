# Last updated march 24 2022

import os
import sys

from Bio import SeqIO

"""This step keeps the bacterial rRNAs that are shared in multiple protists. This
means they are more likely to be symbionts that we can look at more."""


"""Reads the '.cluster' file from the cd-hit-est program. It catches the number
of organisms that have a particular bacterial rRNA sequence."""
def parse_cd_hit_est_clstr(clstr_file):
    raw_data = open(clstr_file).read().split('>Cluster ')[1:]

    cluster_data = {}
    rep_seq_db = {}
    large_format = {}

    # Go through every cluster to keep the representative sequence and all the
    # cells that had the same rRNA.
    for i in raw_data:
        cluster_info = i.split('\n')[:-1]
        taxa_with_seq = []
        cov_info = {}
        cluster_num = cluster_info[0]

        for hit in cluster_info[1:]:
            if '*' in hit:
                rep_seq = hit.split('>')[1].split('...')[0]

            taxon = hit.split('>')[1].split('...')[0]
            taxa_with_seq.append(taxon[:10])

            cov = int(taxon.split('Cov_')[1].split('_')[0])
            cov_info.setdefault(taxon[:10],[]).append(cov)

        # Keep track of how many cells were found with a specific bacterial rRNA.
        unique_taxa_with_seq = [(k, sum(v)/len(v)) for k,v in cov_info.items()]

        rep_seq_db[rep_seq] = f'rRNA_Cluster_{cluster_num}'

        cluster_data[f'rRNA_Cluster_{cluster_num}'] = [rep_seq, unique_taxa_with_seq]

    return cluster_data, rep_seq_db


"""Keeps the bacterial rRNA sequences if they were found in at least 3 different
cells. The minimum number of cells can be changed, 3 was good for the sampling
we had access to."""
def filtering_clusters(cluster_data, min_taxa = 3):

    good_clusters = {}

    for key, value in cluster_data.items():
        number_unique_taxa = len(value[1])

        if number_unique_taxa >= min_taxa:
            good_clusters[key] = value

    return good_clusters


"""Make a table with the name of the cluster, the representative sequence name,
and the number of cells found with the rRNA sequence."""
def write_data_out(cluster_data, good_clusters, rep_seq_db, cluster_fasta, scope = 'Asterionella'):

    # Has all the information for all the clusters.
    cluster_summary = f'{scope}.OverallClusteringSummary.tsv'

    # Has just the final rRNA clusters that were found in at least 3 cells.
    good_cluster_summary = f'{scope}.FinalClusters.Summary.tsv'

    # Table with converted names (to make it easier to keep track).
    rep_seq_conversion = f'{scope}.FinalClustersRenaming.tsv'

    # Sequence file with all the good bacterial rRNAs that are renamed.
    updated_fasta = f'{scope}.Bac_16s_rRNA.Final.fas'

    cluster_seqs = {}
    all_taxa = []
    for seq_rec in SeqIO.parse(cluster_fasta,'fasta'):
        cluster_seqs[seq_rec.id] = f'{seq_rec.seq}'
        all_taxa.append(seq_rec.id[:10])

    all_taxa = list(set(all_taxa))

    with open(cluster_summary,'w+') as w:
        w.write('Cluster_Name\tRepresentative_Sequence\tNumber_of_Cells\n')
        for key, value in cluster_data.items():
            w.write(f'{key}\t{value[0]}\t{len(value[1])}\n')

    with open(good_cluster_summary,'w+') as w:
        w.write('Cluster_Name\tNumber_of_Cells\tRepresentative_Sequence\tTaxon\tCoverage\n')
        for key, value in good_clusters.items():
            rep_seq = value[0]
            total_taxa = len(value[1])
            for i in value[1]:
                w.write(f'{key}\t{total_taxa}\t{rep_seq}\t{i[0]}\t{i[1]}\n')
                # w.write(f'{key}\t{rep_seq}\t{i}\t{total_taxa}\n')

    with open(rep_seq_conversion, 'w+') as w:
        w.write('Original_Name\tNew_Name\n')
        for k, v in rep_seq_db.items():
            if k in good_clusters:
                w.write(f'{k}\t{v}\n')

    with open(updated_fasta, 'w+') as w:
        for k, v in good_clusters.items():
            old_name = v[0]
            new_name = rep_seq_db[old_name]
            rep_seq = cluster_seqs[old_name]
            w.write(f'>{new_name}\n{rep_seq}\n')

if __name__ == '__main__':
    if len(sys.argv[1:]) == 4:
        cluster_file = sys.argv[1]
        cluster_fasta_file = sys.argv[2]
        project = sys.argv[3]
        minimum_cells = int(sys.argv[4])
    else:
        print("\nUsage:\npython filter_clustered_16s.py [Cluster-File] "\
            "[Clustered-FASTA-File] [Project-Name] [Min-Num-Cells]\n")
        sys.exit()

    cluster_data, representative_seqs = parse_cd_hit_est_clstr(cluster_file)
    final_clusters = filtering_clusters(cluster_data, minimum_cells)
    write_data_out(cluster_data,
                final_clusters,
                representative_seqs,
                cluster_fasta_file,
                project)
