from Bio.Blast.NCBIWWW import qblast
from Bio import SeqIO, GenBank
import sys
import re
import pandas as pd
import numpy as np

def get_prots_seqs(name, prots):
    prots_seqs = []
    for record in SeqIO.parse(name, "fasta"):
        prot = record.id
        if prot in prots:
            prots_seqs.append(record)
    return prots_seqs

def load_its(gb_name, species):
    with open(gb_name, 'r') as inp:
        its_species = [x.organism for x in GenBank.parse(inp)]
    its_seqs = list(SeqIO.parse(gb_name,'gb'))
    its_lengths = (len(x) for x in its_seqs)
    its_df = pd.DataFrame({'spec':its_species,'len':its_lengths}).reset_index()
    its_inds = its_df[its_df['spec'].isin(species)].sort_values('len', ascending=False).groupby('spec').first()['index'].values
    its_seqs = [its_seqs[i] for i in its_inds]
    its_species = [its_species[i] for i in its_inds]
    return its_seqs, its_species

def main():
    row = pd.read_excel('results.xlsx', 'BBH').query('puccinia_prots==@PUCCINIA_PROTS')
    puccinia_prots = row['puccinia_prots'].iloc[0].split(',')
    seqs = get_prots_seqs('puccinia_graminis.fasta', puccinia_prots)
    neurospora_prots = row['neurospora_prots'].iloc[0].split(',')
    seqs = seqs + get_prots_seqs('neurospora_crassa.fasta', neurospora_prots)
    
    others_seqs = list(SeqIO.parse(PUCCINIA_PROTS + '.txt', "fasta"))
    species = [re.findall('\[(.*)\]',x.description)[0] for x in others_seqs]
    others_seqs, species = list(map(list, zip(*sorted(zip(others_seqs, species), key=lambda x: x[1]))))
    _, index = np.unique(species, return_index=True)
    others_seqs = others_seqs[:index[23]]
    seqs = seqs + others_seqs
    species = [' '.join(re.findall('\[(.*)\]',x.description)[0].split(' ')[:2]) for x in seqs]
    query = '(' + ' OR '.join(['"' + s.replace('"','-') + '"[Organism]' for s in set(species)]) + ')' + 'AND "internal transcribed spacer"[TI]'
    print(query) # use that to get its.gb from ncbi
    its_seqs, its_species = load_its('its.gb', species)


    for seq, spec in zip(seqs, species):
        seq.id =  '_'.join([x.strip() for x in spec.split(' ')]) + "_" +  seq.id
    for seq, spec in zip(its_seqs, its_species):
        seq.id =  '_'.join([x.strip() for x in spec.split(' ')])
    with open('phylo_prot.fasta', 'w') as out:
        SeqIO.write(seqs, out,'fasta')
    with open('phylo_its.fasta', 'w') as out:
        SeqIO.write(its_seqs, out,'fasta')  
    pass

PUCCINIA_PROTS = 'XP_003321789.2'
if __name__ == '__main__':
    if len(sys.argv) > 1:
        PUCCINIA_PROTS = (sys.argv[1])
    main()
