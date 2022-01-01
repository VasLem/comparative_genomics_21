from Bio import SeqIO
from Bio.SearchIO.BlastIO.blast_xml import BlastXmlParser
from Bio.SeqIO import parse
from Bio import Entrez
import pandas as pd
import re
from urllib.request import urlretrieve 
from tqdm.auto import tqdm
tqdm.pandas()
from tqdm import tqdm
import os


Entrez.email = 'vasileios.lemonidis@student.kuleuven.be'
def download_sequence(row, minus_tss=3000, plus_tss=3000):
    os.makedirs('re_regions', exist_ok=True)
    fname = os.path.join('re_regions', f"{row['query'].split()[0]}_{minus_tss}_{plus_tss}.fasta")
    if not os.path.isfile(fname):
        strand = 2 if row['strand'].lower() == "minus" else 1
        handle = Entrez.efetch(db="nucleotide", rettype="fasta", id=row['gi_id'],
                               seq_start=row['start_pos'] + 1 - minus_tss,
                               seq_stop=row['start_pos'] + 1 + plus_tss,
                               strand=strand)                       
        ret = SeqIO.read(handle, 'fasta')
        ret.id = '_'.join(row['organism'].split()) + '_' + row['query'].split()[0] + '_' + ret.id
        ret.description = row['query'].split()[0]
        SeqIO.write(ret, fname, 'fasta')
        return ret
    else:
        return SeqIO.read(fname, 'fasta')

    

# import gzip
# import shutil
# def find_sequence(row):
#     genome = next(parse(os.path.join('whole_genomes',row['organism']),'fasta'))
#     if row['strand'] == 'minus':
#         genome.seq = genome.seq[::-1].complement()
#     return genome.seq[row['start_pos'] - 100: row['start_pos'] + 50]
# def download_genomes(organisms):
#     eukaryotes_genomes = pd.read_csv('eukaryotes.csv', sep=',')#.dropna(subset=["RefSeq FTP"])
#     mask = [any(y.startswith(x) for y in organisms) for x in eukaryotes_genomes['#Organism Name']]
#     eukaryotes_genomes = eukaryotes_genomes.loc[mask,:]
#     eukaryotes_genomes['query'] = [[y for y in organisms if y.startswith(x)][0] for x in eukaryotes_genomes['#Organism Name']]
#     os.makedirs('whole_genomes', exist_ok=True)
#     for _, row in tqdm(eukaryotes_genomes.iterrows()):
#         fpath = os.path.join('whole_genomes', row['query'])
#         if not os.path.exists(fpath):
#             catalog_fil = urlretrieve(row['RefSeq FTP'])[0]
#             with open(catalog_fil,'r') as inp:
#                 catalog = inp.read()
#             fname = re.findall(r"(\S*\_genomic\.fna\.gz)", catalog)
#             fname = [x for x in fname if '_rna_' not in x and '_cds_' not in x][0]
#             gz_fil = urlretrieve(os.path.join(row['RefSeq FTP'], fname))[0]
#             with gzip.open(gz_fil, 'rb') as f_in:
#                 with open(fpath, 'wb') as f_out:
#                     shutil.copyfileobj(f_in, f_out)

def download_loci_info(genes_ids):
    handle = Entrez.efetch('gene',id=','.join(genes_ids), retmode="xml")
    records = list(Entrez.parse(handle))
    starts = []
    ends = []
    strands = []
    gids = []
    for r in records:
        loc_r =  [x for x in r['Entrezgene_locus'] if 'Gene-commentary_seqs' in x][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']
        starts.append(int(loc_r['Seq-interval_from']))
        ends.append(int(loc_r['Seq-interval_to']))
        gids.append(loc_r['Seq-interval_id']['Seq-id']['Seq-id_gi'])
        strands.append(loc_r['Seq-interval_strand']['Na-strand'].attributes['value'])
    return starts, ends, gids, strands

def download_genes_ids(rna_accessions):
    handle = Entrez.efetch('nuccore',id=','.join(rna_accessions), retmode="xml")
    records = list(Entrez.parse(handle))
    genes = []
    for r in records:
        gene_r = [x for x in r['GBSeq_feature-table'] if x['GBFeature_key']=='gene'][0]['GBFeature_quals']
        gene_id = [g for g in  gene_r if g['GBQualifier_name'] == 'db_xref'][0]['GBQualifier_value'].split(":")[1]
        genes.append(gene_id)
    return genes


def main():
    accessions = []
    queries = []
    blast_scores = []
    with open('prot2rna.xml','r') as inp:
        for t in BlastXmlParser(inp):
            first_hit = t[0]
            accessions.append(first_hit.accession)
            queries.append(first_hit.query_description)
    organisms = [re.findall(r'\[(.*)\]', q)[0] for q in queries]
    # download_genomes(organisms)
    genes = download_genes_ids(accessions)
    starts, ends, gi_ids, strands = download_loci_info(genes)
    genetic_info = pd.DataFrame({'query':queries, 'organism':organisms, 'start_pos':starts, 'end_pos':ends, 'gi_id': gi_ids, 'strand':strands})
    genetic_info['sequence'] = genetic_info.progress_apply(download_sequence, axis=1)
    SeqIO.write(genetic_info['sequence'].to_list(), 'regulatory_regions.fasta', 'fasta')
    with pd.ExcelWriter("results.xlsx",mode='a') as writer:
        workBook = writer.book
        try:
            workBook.remove(workBook['loci'])
        except:
            pass
        finally:
            genetic_info.to_excel(writer, sheet_name='loci', index=False)
            writer.save()

    


if __name__ == '__main__':
    main()
