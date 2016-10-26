import sarge
import re
import os

def Honey_pie(sortBam,sortTailBam,ref_fa,thread,tmp,otherParams=['']):
    """Honey pip extract soft clip reads and remap them
    """
    tailBam = re.sub('\.final\.bam$','.tail.bam',sortTailBam)
    cmd = ('Honey.py pie -o {tail} -n {thread} {input} {ref} --temp {tmp}').format(tail=tailBam,
                            thread=str(thread),input=sortBam,ref=ref_fa,tmp=tmp)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)
    # sort
    cmd = ('samtools sort -m 4G -@ {thread} -T {pre} -o {sortBam} {bam} ').format(
            thread=str(thread),pre=tailBam[:-4],sortBam=sortTailBam,bam=tailBam)
    print(cmd)
    sarge.run(cmd)
    # index
    cmd = ('samtools index {out} ').format(out=sortTailBam)
    print(cmd)
    sarge.run(cmd)
#     os.remove(sortBam)
    

def Honey_tails(finalBam,bamTail,otherParams=['']):
    """This function run Honey tail,culster the soft clipped reads
    """
    cmd = ('Honey.py tails -o {out} {input} ').format(input=finalBam,out=bamTail)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)


def Honey_spots(finalBam,spotFile,ref_fa,thread,otherParams=['']):
    """This function run Honey sorts.
    """
    cmd = ('Honey.py spots --reference {ref} -n {thread} -o {out} {input} ').format(
            input=finalBam,ref=ref_fa,thread=str(thread),out=spotFile)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)

import pandas as pd
import numpy as np
class pb_tail_res(object):
    '''
    Input is output result from pb tail. should be pandas dataframe
    '''
    def __init__(self,df):
        self.df = df
        self.df.columns = ['id','chrKey','uRef','uBreak','uMapq','dRef','dBreak','dMapq','remainSeq','annot','numReads','numZMWs','evidence']
    
    def get_sv_types(self):
        '''get all the sv types'''
        types = list(set(self.df['annot'].tolist()))
        return types
        
    def get_sv_num(self,sv_type):
        '''get sv number'''
        df = self.df
        return df[df['annot'].values==sv_type].shape[0]
    
    def add_sv_len(self):
        '''add sv length for each sv except translocation whose break points are in different chromosome.'''
        df = self.df
        df['len'] = df.apply(lambda row: 'NA' if row['annot']=='TLOC' else int(row['dBreak'])-int(row['uBreak']),axis=1)
        return df
    
    def get_sv_num4_each_chr(self,chr_len_df,sv_type,count_log=False,length_log=False):
        '''get sv number for each scaffold for specific sv type
        * chr_len_df: pandas dataframe with 1 column. ['chr_len']. chr name is index
        '''
        df = self.df[self.df['annot'].values==sv_type]
        sv_count = df.groupby(['uRef']).size()
        df = pd.concat([sv_count,chr_len_df],axis=1)
        df = df.fillna(0)
        df = df.rename(columns={0:'count'})
        if count_log == True:
            df['count'] = np.log10(df['count'])
        if length_log == True:
            df['chr_len'] = np.log10(df['chr_len'])
        return df


    
        
        
        