import re
import pandas as pd


class ncbi_gff(object):
    def __init__(self,df):
        self.df = df
        self.df.columns=['chr','source','feature','start','end','score','strand','frame','anno']
        self.df['start'] = self.df['start'] - 1
        self.df = self.df[self.df['feature'].values!='region']
        self.df = self.df.reset_index(drop=True)
        self.df['geneid'] = self.df['anno'].apply(lambda x: ncbi_gff.get_id(x,'GeneID:'))
        self.df['trid'] = self.df['anno'].apply(lambda x: ncbi_gff.get_id(x,'transcript_id='))
        self.df['prid'] = self.df['anno'].apply(lambda x: ncbi_gff.get_id(x,'protein_id='))
    @staticmethod
    def get_id(anno,feature):
        '''get id based on the feature provided'''
        try:
            gene_id = re.search('(?<={id}).+?(?=[;,]|$)'.format(id=feature),anno).group(0)
        except:
            gene_id = None
        return gene_id
    
    @staticmethod
    def get_tr_longest_intron(tr_df):
        '''get the longest intron the the transcript'''
        start = tr_df['start'].tolist()
        end = tr_df['end'].tolist()
        strand = tr_df['strand'].tolist()
        if len(start) == 1:
            return 0
        if strand[0] == '+':
            intron = max([abs(int(s)-int(e)) for s,e in zip(start[1:],end[:-1])])
        else:
            intron = max([abs(int(s)-int(e)) for s,e in zip(start[:-1],end[1:])])
        return intron
    
    def get_longest_intron(self):
        '''this is the longest intron across the whole genome'''
        df = self.df
        df = df[(~df['prid'].isnull()) & (df['feature'].values=='CDS')]
        df = df.reset_index(drop=True)
        df = df.groupby(['chr','prid']).apply(ncbi_gff.get_tr_longest_intron)
        return df#.max()
    
    def get_all_id(self):
        '''this function gets all ids in the gff file
        '''
        df = self.df
        id_df = df[df['feature'].isin(['exon','CDS'])]
        id_df = id_df.reset_index(drop=True)
        
        id_df['sym'] = id_df['anno'].map(lambda x: ncbi_gff.get_id(x,'gene='))
        id_df['rna'] = id_df['anno'].map(lambda x: ncbi_gff.get_id(x,'Parent='))
        
        exn_df = id_df[id_df['feature'].values=='exon'][['geneid','sym','chr','rna','trid']].drop_duplicates()
        exn_df = exn_df.reset_index(drop=True)
        
        cds_df = id_df[id_df['feature'].values=='CDS'][['geneid','sym','chr','rna','prid']].drop_duplicates()
        cds_df = cds_df.reset_index(drop=True)
        
        merge_df = pd.merge(exn_df,cds_df,how='outer',on=['geneid','sym','chr','rna'])
        merge_df.columns = ['GeneID','GeneSymbol','Chrom','TrID','TrAccess','PrAccess']
        merge_df.fillna('-',inplace=True)
        merge_df.fillna('-',inplace=True)
#         merge_df = merge_df[(merge_df['TrAccess'].values != '-') | (merge_df['PrAccess'].values != '-')]
        merge_df = merge_df.sort_values(['GeneID'])
        return merge_df[['GeneID','GeneSymbol','Chrom','TrAccess','PrAccess','TrID']]
    
    def get_gene_seq(self,ref_dic,gid,id_type='tr'):
        '''this function gets seqeunce of a transcript or protein
        '''
        df = self.df
        if id_type == 'tr':
            feature = 'exon'
            id_t = 'trid'
        elif id_type == 'pr':
            feature = 'CDS'
            id_t = 'prid'
        region_df = df[(df['feature'].values==feature) & (df[id_t].values==gid)]
        # get sequence
        scaff = region_df['chr'].tolist()[0]
        scaff_seq = ref_dic[scaff].seq
        strand = region_df['strand'].tolist()[0]
        
        g_seq = ''
        for s,e in zip(region_df['start'],region_df['end']):
            g_seq += scaff_seq[int(s):int(e)]
        # consider strand
        if strand == '-':
            g_seq = g_seq.reverse_complement()
        
        if id_type == 'pr':
            g_seq = g_seq.translate()
        return g_seq
        
# gff_fn = '/data/genome/hamster/ncbi_refseq/hamster.gff'
# df = pd.read_csv(gff_fn,sep='\t',header=None,comment='#')
# obj = ncbi_gff(df)
# all_id_df = obj.get_all_id()
# all_id_df.to_csv('/data/genome/hamster/ncbi_refseq/all_id.txt',sep='\t',index=False)
# 
# ref_dic = SeqIO.index('/data/genome/hamster/picr/picr.fa','fasta')
# res = obj.get_gene_seq(ref_dic,'NM_001246795',id_type='tr')
# print res


from Bio import SeqIO

def getr_tRNA(fa,gff,output):
    '''get rRNA and tRNA sequence
    * output: fa file stores rtRNA sequence'''
    index = SeqIO.index(fa,'fasta')
    with open(gff) as f, open(output,'w') as out:
        for line in f:
            if line.startswith('#'):continue
            item = line.strip().split('\t')
            if item[2] in ['rRNA','tRNA']:
                name = re.search('(?<=product=).+?(?=$|;)',line).group(0)
                chrom = item[0]
                s = int(item[3])
                e = int(item[4])
                seq = str(index[chrom].seq[s-1:e])
                out.write('>'+name + '\n' + seq + '\n')