# Last updated April 6, 2022

"""Cleans up the table outpus from the Vsearch tool. It keeps just the deeper
taxonomy (Phylum, Class, Order, Family) since the data look bad at the species level."""

import os
import sys
import pandas as pd


"""Keep just the useful parts of the table and flip it so each bacteria has a row,
while the expression information for each cell are in the columns."""
def clean_table(rRNA_clust_table, rRNA_assign_table, project_name, num_cells = 3):
    # rRNA cluster table that has the data for each protist cell.
    rRNA_clust = pd.read_table(rRNA_clust_table)

    # Table with the taxonomy information from VSearch.
    rRNA_assign = pd.read_table(rRNA_assign_table)

    # Keep just the information of the good clusters that we could assign to
    # a known bacterium.
    clean_table = rRNA_clust[rRNA_clust["Cluster_Name"].isin(rRNA_assign["Cluster_Name"])]

    # Setting up our table for the next steps.
    fin_table = clean_table.pivot_table(index = 'Cluster_Name',
                                        columns = 'Taxon',
                                        values = 'Coverage')

    taxonomy = {row['Cluster_Name']:get_taxonomy(row['rRNA_Hit']) for index, row in rRNA_assign.iterrows()}

    # Keep just the basic taxonomy information as we don't know much about strains.
    temp_df = pd.DataFrame.from_dict(taxonomy, columns = ['Phylum','Class','Order','Family'], orient='index')
    temp_df.index.name = 'Cluster_Name'

    # Fill missing data as "0" to make visualizing the data easier.
    fin_df = fin_table.fillna(0).astype('int')

    finished_df = temp_df.merge(fin_df, how='inner', on='Cluster_Name')

    # Save the final table!
    finished_df.to_csv(f'{project_name}.Microbiome.{num_cells}Cells.BigTable.tsv', sep = '\t')


"""Breaks down the taxonomy from Vsearch and saves just the useful parts."""
def get_taxonomy(rRNA_hit):

    hit_name = rRNA_hit.split('_')[1].split(';')[1:]

    t_phylum = hit_name[0]
    t_class = hit_name[1]
    t_order = hit_name[2]
    t_family = hit_name[3]

    return [t_phylum, t_class, t_order, t_family]


if __name__ == '__main__':
    if len(sys.argv[1:]) == 3:
        rRNA_clust_table = sys.argv[1]
        rRNA_assign_table = sys.argv[2]
        project_name = sys.argv[3]

    else:
        print("\nUsage:\npython vsearch_16s.py [cluster-table] [assignment-table] "
            "[PROJECT NAME]\n")
        sys.exit()

    clean_table(rRNA_clust_table, rRNA_assign_table)
