# Last updated march 31 2022

"""Uses the vsearch tool to search for the taxonomy of our rRNA sequences."""

import os
import sys

"""Runs the vsearch tool against the database (SilvaDB)."""
def run_vsearch(query, database, project_name, identity=80):
    vsearch_cmd = f'vsearch --usearch_global {query} ' \
        f'--db {database} ' \
        f'--id {identity/100} ' \
        f'--maxaccepts 1 ' \
        f'--blast6out {project_name}.16s_Assessment.{identity}id.tsv'

    os.system(vsearch_cmd)

    clean_table(project_name, identity)

"""Cleans the outputs from the vsearch tool, and adds a header line to make it
easier to read the table."""
def clean_table(project_name, identity=95):

    info = [i for i in open(f'{project_name}.16s_Assessment.{identity}id.tsv').readlines()]

    with open(f'{project_name}.16s_Assessment.{identity}id.tsv', 'w+') as w:
        w.write('Cluster_Name\trRNA_Hit\tIdentity\tAln_Length\tMismatches\tGaps\tQ_St\tQ_End\tH_St\tH_End\tEval\tHSP\n')
        w.write(''.join(info))



if __name__ == '__main__':
    if len(sys.argv[1:]) == 3:
        query_16s_fasta = sys.argv[1]
        database_fasta = sys.argv[2]
        project_name = sys.argv[3]
    else:
        print("\nUsage:\npython vsearch_16s.py [QUERY FASTA] [DATABASE FASTA] "
            "[PROJECT NAME]\n")
        sys.exit()

    run_vsearch(query_16s_fasta, database_fasta, project_name, identity=95)
