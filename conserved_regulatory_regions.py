import os
import pandas as pd
from pyjaspar import jaspardb
from uniprot import batch_uniprot_metadata

from conserved_regions import analyze_conserved
from utils import append_sheet_to_excel
import requests


def get_uniprot_id(row):
    jdb_obj = jaspardb(release='JASPAR2022')
    if row['Jaspar motif hit'].strip() == '':
        return None
    motif = jdb_obj.fetch_motif_by_id(row['Jaspar motif hit'].split()[0])
    return motif.acc[0]

from io import StringIO
import time

from bioservices import UniProt
import re

def get_uniprot_function(row):
    u = UniProt(verbose=False)
    res = u.search(row['Uniprot ID'], frmt="tab", limit=3, columns="comment(FUNCTION)")
    #res = u.retrieve(value, frmt='xml', database='uniparc')
    return re.findall('FUNCTION:(.*)', res)[0]

def main():
    analyze_conserved('regulatory_regions.aln', consensus_threshold=0.7,smooth_size=11, conserve_threshold=0.45, conserved_size_min=6)
    if not os.path.isdir("meme_output"):
        os.system('meme-chip regulatory_regions.fasta -db JASPAR2022_CORE_fungi_non-redundant_pfms_meme.txt  -minw 8 -maxw 30  -oc meme_output -dna -meme-nmotifs 20 --ccut 0')
    summary = pd.read_csv('meme_output/summary.tsv', sep='\t',comment='#')[['CONSENSUS', 'E-VALUE', 'E-VALUE_SOURCE', 'MOST_SIMILAR_MOTIF']]
    summary.columns = ['Motif', 'Hit E-value', 'MEME tool', 'Jaspar motif hit']
    summary['Uniprot ID'] = summary.apply(get_uniprot_id, axis=1)
    desc_df = summary.dropna(subset=['Uniprot ID'])[['Motif', 'Uniprot ID']]
    desc_df['Uniprot function'] = desc_df.apply(get_uniprot_function, axis=1)
    summary = summary.merge(desc_df, on='Motif', how='left')
    append_sheet_to_excel('results.xlsx', 'promoter_regions_motifs', summary, index=False)
    

    

if __name__ == '__main__':
    main()