from Bio.AlignIO import read
from Bio.Align.AlignInfo import SummaryInfo
from Bio.SeqUtils import IUPACData
import numpy as np
import matplotlib.pyplot as plt
from seaborn import rugplot
from numpy.lib.stride_tricks import sliding_window_view
from matplotlib.collections import LineCollection
import pandas as pd
import os
from io import StringIO
from utils import append_sheet_to_excel
from batch_cd_search import search
import re
def analyze_conserved(name, consensus_threshold=0.7, conserve_threshold=0.3, conserved_size_min=8, criterion='Simpson diversity', smooth_size=11):
    """
    consensus_threshold: minimum frequency for aa to be considered of high certainty
    conserve_threshold: maximum criterion threshold for aa to be considered conserved
    conserved_size_min: minimum size of a conserved sequence
    criterion: "Shannon entropy" or "Simpson diversity"[default]
    smooth_size: smoothening factor for amino acids sequences, needs to be odd number or 0, only applicable to the plot
    """
    assert criterion in ["Shannon entropy", "Simpson diversity"]
    assert smooth_size % 2 == 1
    aas = {a:c + 1 for c,a in enumerate(IUPACData.protein_letters)}
    aas['-'] = 0
    aligned = read(name, 'clustal')
    aligned_info = SummaryInfo(aligned)
    consensus = aligned_info.gap_consensus(threshold=consensus_threshold)
    seqs = list(map(list, zip(*[[ aas[y] for x in str(t.seq) for y in x] for t in aligned])))
    counts = [np.unique(x, return_counts=True,)[1] for x in seqs]
    probs = [x/np.sum(x) for x in counts]
    if criterion == 'shannon':
        criterion_val = np.array([-np.sum(x * np.log2(x)) for x in probs])
        criterion_val = criterion_val / np.max(criterion_val)
    else:
        criterion_val = np.array([1 - np.sum(x*(x-1))/(np.sum(x) * (np.sum(x) - 1)) for x in counts])
    
    smoothen = lambda x: np.mean(sliding_window_view(np.pad(x, [smooth_size//2,smooth_size//2], 'edge'), smooth_size),axis=1) if smooth_size > 0 else x
    fig, ax = plt.subplots(figsize=(15,15))
    xy = np.array([np.arange(len(criterion_val)), smoothen(criterion_val)]).T
    xy = xy.reshape(-1, 1, 2)
    segments = np.hstack([xy[:-1], xy[1:]])
    missing = [np.mean(np.array(s)==0) for s in seqs]
    coll = LineCollection(segments, cmap=plt.cm.gnuplot, label=f'{criterion}')
    coll.set_array(np.array(missing))
    ax.add_collection(coll)
    ax.autoscale_view()
    
    mask = criterion_val < conserve_threshold
    
    axcb = fig.colorbar(coll)
    axcb.set_label('Gaps Ratio')
    ax.set_title(f"Sequence Conservation")
    ax.hlines(conserve_threshold,xmin=0,xmax=len(criterion_val),label='Threshold', colors=['red'])
    ax.set_xlim([0, len(mask)])
    plt.legend()
    diff = np.diff(np.pad(mask.astype(int), [1,1], 'constant'))
    # conserved regions of minimum size conserved_size_min and not containing uncertain concensus or gaps
    conserved_regions =  [x for x in list(map(list, zip(*[np.where(diff==1)[0] , np.where(diff==-1)[0] - 1]))) if
                          x[1]-x[0] >=  conserved_size_min and
                          '-' not in consensus[x[0]:x[1]+1].replace('X','-')]
    exp_mask = np.zeros(len(mask))
    for x in conserved_regions:
        exp_mask[x[0]:x[1]+1] = 1
    
    if conserved_regions:
        result = pd.DataFrame(conserved_regions)
        result.columns = ['start', 'end']
        result['sequence'] = [''.join(consensus[s:e+1]) for s,e in conserved_regions]
        cd_search = search([str(x) for x in result['sequence'].to_list()])
        cd_search = cd_search[[c for c,x in enumerate(cd_search) if x.startswith('Query')][0]:]
        cd_search = '\n'.join(cd_search).rstrip()
        cd_ret = pd.read_csv(StringIO(cd_search),sep='\t')
        cd_ret = cd_ret[cd_ret['Hit type'] == 'superfamily'].drop(columns=['Hit type', 'Bitscore', 'Superfamily'])
        cd_ret.drop(columns=['From', 'To'],inplace=True)
        cd_ret = cd_ret.sort_values('E-Value').groupby("Query").first().reset_index()
        cd_ret['sequence'] = cd_ret['Query'].apply(lambda x: result['sequence'][int(re.findall('>(.*)',x)[0])])
        cd_ret.drop(columns='Query',inplace=True)
        cd_ret.rename(columns={'Short name':'Superfamily', 'Accession':'Superfamily accession'},inplace=True)
        cd_ret['Superfamily'] = cd_ret['Superfamily'].apply(lambda x: x.replace('superfamily', ''))
        result = result.merge(cd_ret, how='left', on='sequence')
        append_sheet_to_excel('results.xlsx', f'{os.path.splitext(name)[0]}conserved', result, index=False)

        
        rugplot(x=np.where(exp_mask)[0])
    
    plt.show()
    fig.savefig(name.split(os.sep)[0] + '.png')
    pass
if __name__ == '__main__':
    analyze_conserved('phylo_prot.aln')
    