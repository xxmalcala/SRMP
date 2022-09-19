# Last updated April 7, 2022

"""Helpful plotting tools to visualize our protist microbiome data. These
functions do not work on their own and are intended for interactive work
in the Python interpreter (in the terminal)."""

import seaborn as sns
import pandas as pd
from venn import venn
import matplotlib.pyplot as plt


"""This function will not work on its own. You must update the lines/locations
and run interactively in the python interpreter."""
def make_venn_diagram():
    # these are the steps to make the venn diagram!
    tsv_table = 'Path/to/BigTable.VennD.tsv'

    tsv_lines = [line.rstrip() for line in open(tsv_table).readlines()]

    rRNA = {'Sr_st_Aformosa':[], 'Sr_st_Ptricorn':[]}

    for line in tsv_lines[1:]:
        cluster = line.split('\t')[0]
        hele = line.split('\t')[1]
        hpap = line.split('\t')[2]
        if hele != '0':
            rRNA['Sr_st_Aformosa'].append(cluster)
        if hpap != '0':
            rRNA['Sr_st_Ptricorn'].append(cluster)

    for k, v in rRNA.items():
        rRNA[k] = set(v)

    venn(rRNA)
    plt.show()
    # Close the window of the plot to get back to the interactive interpreter!


"""Makes a heatmap using the "BigTable.tsv" from step 6_finalize_table.py"""
def make_heatmap():
    # heatmaps with Seaborn: https://seaborn.pydata.org/generated/seaborn.heatmap.html
    tsv_table = 'Path/to/BigTable.tsv'
    df = pd.read_table(tsv_table)

    # Removing columns from the table that we don't want to work with
    df2 = df.drop(columns=['Class', 'Phylum', 'Order', 'Family'])

    # Set the y-axis here! Whatever is in this column will become the axis
    # labes.
    df3 = df2.set_index('Cluster_Name')

    # Some snippet to help with masking parts of the heatmap based on
    # their numerical values
    values = df3.to_numpy(dtype=float)

# Heatmap command! Cmap --> Color map, vmin/vmax --> Min/Max values for the scaling,
# mask = values > number --> Mask parts of the heatmap below some value!
ax = sns.heatmap(df3, xticklabels=True, cmap = 'coolwarm', yticklabels=True, vmin=10, vmax=200, mask=values < 10, linewidths=0.5, linecolor='lightgrey')
for _, spine in ax.spines.items():
    spine.set_visible(True)

# ax = sns.clustermap(df3, xticklabels=True, cmap = 'coolwarm', yticklabels=True, vmin=10, vmax=200, mask=values < 10, linewidths=0.5, linecolor='lightgrey')
