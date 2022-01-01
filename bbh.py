import pandas as pd
import re
from Bio import SeqIO
import csv
import os
from utils import append_sheet_to_excel

def load_df(name):
    fields_str = r'query, target, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, identity'
    cols = fields_str.split(', ')
    df = pd.read_csv(name, sep='\t', header=None,index_col=None)
    df.columns = cols
    df = df[['query', 'target', 'identity']].copy()
    return df

def preprocess_self_blast(df: pd.DataFrame):
    #rename the columns
    df.rename(columns={'identity':"self_identity", 'target':'self_target'}, inplace=True)
    #remove identity hits 
    df = df[df['self_target'] != df['query']]
    # Keep only the query from the strain with the best hit
    df = df.sort_values('self_identity',ascending=False).groupby('query').first().reset_index()
    return df

def preprocess_cross_blast(df: pd.DataFrame):
    # Keep only the query from the strain with the best hit
    df = df.sort_values('identity',ascending=False).groupby('query').first().reset_index()
    return df


def find_inparalogs(self_df, df):
    ret = self_df[['query', 'self_identity', 'self_target']].merge(df, on='query')
    # in-paralog : when the best hit score inside the same species is less than the best hit score with the different species
    ret['is_inparalog'] = (ret['self_identity'] != 100)  & (ret['self_identity'] > ret['identity'])
    # group all in-paralogs that share the same best hit in-species target
    inparalogs = ret[ret['is_inparalog']].groupby('self_target')['query'].aggregate(list)
    inparalogs.name = 'inparalogs'
    inparalogs.index.name = 'query'
    df = df.merge(inparalogs, on='query', how='left')
    df['inparalogs'] = df['inparalogs'].fillna(df['query'].apply(lambda x: []))
    # add to the inparalogs list the query itself (which was the target for the in-paralogs)
    df['inparalogs'] = df.apply(lambda x: tuple(set([x['query']] + x['inparalogs'])), axis=1)
    # drop duplicates by target and in-paralogs, and keep the one with the highest recorded score
    df = df.sort_values('identity',ascending=False).groupby(['target', 'inparalogs']).first().reset_index()
    return df

def main():
    q_puccinia = preprocess_cross_blast(load_df('query_puccinia_graminis.out'))
    s_puccinia = preprocess_self_blast(load_df('self_puccinia_graminis.out'))
    puccinia_df = find_inparalogs(s_puccinia, q_puccinia)

    q_neurospora = preprocess_cross_blast(load_df('query_neurospora_crassa.out'))
    s_neurospora = preprocess_self_blast(load_df('self_neurospora_crassa.out'))
    neurospora_df = find_inparalogs(s_neurospora, q_neurospora)

    merged_df = puccinia_df.merge(neurospora_df, left_on='query', right_on='target', suffixes=['_puc','_neu'], how='inner')
    
    bbh = merged_df.query('target_puc==query_neu')
    bbh = bbh[['inparalogs_puc', 'inparalogs_neu', 'identity_puc']].rename(
        columns={'inparalogs_puc':"puccinia_prots",
                 'inparalogs_neu':"neurospora_prots",
                 "identity_puc": 'identity'}).groupby(
                    ['puccinia_prots','neurospora_prots']).aggregate(max).reset_index()
    report = pd.DataFrame()

    report.loc['puccinia graminis','GenesNumber'] = len(set(q_puccinia['query']))
    report.loc['neurospora crassa','GenesNumber'] = len(set(q_neurospora['query']))
    report.loc['puccinia graminis','OrthologousNumber'] = len(set([x for y in bbh['puccinia_prots'] for x in y]))
    report.loc['neurospora crassa','OrthologousNumber'] = len(set([x for y in bbh['neurospora_prots'] for x in y]))

    bbh['puccinia_len'] = bbh['puccinia_prots'].apply(len)
    bbh['neurospora_len'] = bbh['neurospora_prots'].apply(len)
    bbh['puccinia_prots'] = bbh['puccinia_prots'].apply(lambda x: ','.join(x))
    bbh['neurospora_prots'] = bbh['neurospora_prots'].apply(lambda x: ','.join(x))

    append_sheet_to_excel('results.xlsx', 'orthologs', bbh, index=False)
    append_sheet_to_excel('results.xlsx', 'orthologs_report', report, index=True)
    
if __name__ == '__main__':
    main()