# coding: utf8
import sarge,sys,os
import pandas as pd 
import numpy as np
import matplotlib as mpl 
mpl.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#Functions to compare results of start-site sequencing to annotation file
#Uses peaks tsv files generated from Homer along with the genome annotation gff3 and gtf files

def add_fields_new_annotation(annotation):
    genome_anno = pd.read_csv(annotation,sep='\t',index_col=0)
    genome_anno['gene'] = genome_anno.index.str.extract('gene=(.*);transcript_id*',expand=False)
    genome_anno['gene_id'] = genome_anno.index.str.extract('gene_id=(.*);gene*',expand=False)
    genome_anno['transcript_id'] = genome_anno.index.str.extract('gene_id=(.*);gene*',expand=False)
    return genome_anno

    
def peaks_within_x(peaks_file,distance=1000):
    curr_peaks = pd.read_csv(peaks_file,sep='\t',comment='#')
    curr_peaks.dropna(axis=1,how='all',inplace=True)
    
    curr_peaks.dropna(axis=0,inplace=True)
    curr_peaks = curr_peaks[np.abs(curr_peaks['Distance to TSS']) < distance]
    curr_peaks['gene'] = curr_peaks['Nearest PromoterID'].str.extract('gene=(.*);transcript_id*',expand=True)
    curr_peaks['gene_id'] = curr_peaks['Nearest PromoterID'].str.extract('gene_id=(.*);gene*',expand=True)
    curr_peaks['transcript_id'] = curr_peaks['Nearest PromoterID'].str.extract('gene_id=(.*);gene*',expand=True)
    return curr_peaks

def gene_centric_TSS(curr_peaks,annotation,divergent_txn_len = 200):
    gene_TSS = pd.DataFrame(index=set(curr_peaks['gene_id'].dropna()),
             columns=['hasPromoter','hasExon1','hasIntron1','XhasExon1','XhasIntron1','Number Of Promoters',
                     'Divergent Transcripts','Introns','Exons','Left Most TSS',
                     'TSS Peak Indices','Intergenic Peaks Adjacent of Gene','has Intergenic Peaks Adjacent of Gene'])
    
    
    gene_TSS.index.names = ['ID']
    gene_TSS = pd.concat((annotation,gene_TSS),axis=1,join='inner')

    
    genes_grouped = curr_peaks.groupby('gene_id')

    counter = 0
    for ind,val in tqdm.tqdm(gene_TSS.iterrows()):
        counter += 1
    #    print(ind)
        tmp = genes_grouped.get_group(ind)
        tmp_promoters = tmp  #tmp[tmp['Reduced Annotation'] == 'promoter']
        tmp_pos_promoters = tmp_promoters[tmp_promoters['Strand'] == val['Strand']]
        tmp_neg_promoters = tmp_promoters[tmp_promoters['Strand'] != val['Strand']]

        if val['Strand'] == '-':
            strand_norm = -1
        elif val['Strand'] == '+':
            strand_norm = 1
        else:
            print('Strand not proper!')

        num_divergent_transcripts = 0
        for ind2,val2 in (tmp_pos_promoters.iterrows()):
            if np.sum((0 < strand_norm*(val2['Start'] - tmp_neg_promoters['End'])) 
                        & (strand_norm*(val2['Start'] - tmp_neg_promoters['End']) <= divergent_txn_len)) > 0:
                num_divergent_transcripts += 1
        num_SS = np.sum(tmp_promoters['Strand'] == val['Strand'] )
        gene_TSS.set_value(ind,'hasPromoter',num_SS > 0)
        gene_TSS.set_value(ind,'Number Of Promoters',num_SS)
        gene_TSS.set_value(ind,'Divergent Transcripts',num_divergent_transcripts)
        gene_TSS.set_value(ind,'hasDivergent',num_divergent_transcripts > 0)        
        has_exon1 = np.sum(tmp['Annotation'].str.contains('exon 1')) > 0
        has_intron1 = np.sum(tmp['Annotation'].str.contains('intron 1')) > 0
        
    return gene_TSS