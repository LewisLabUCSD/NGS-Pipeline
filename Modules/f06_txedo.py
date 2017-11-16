import subprocess,os
from natsort import natsorted
import pandas as pd
"""
this file includes all functions about txedo pipeline. (cufflinks,cufmerge,cuffdiff) 
"""
def cufflinks(bamFiles,annotation,thread,otherParameters=['']):
    """
    this function run cufflinks.
    
    * bamFiles:           a list of bam files. [1.bam,2.bam,...]
    * annotation:         gtf or gff files
    * thread:             umber of cores
    * otherParameters:    additional parameters in a list [param1,param2,...]
    """
    cmd = ''
    result = []
    for bam in bamFiles:
        out_dir = bam.split('.')[0] +'_cufflinks'
        result.append(out_dir)
        cufflinkCmd = ('cufflinks -o {out_dir} -p {thread} -G {annotation} {bam}').format(
                        out_dir=out_dir,thread=thread,annotation=annotation,bam=bam)
        cmd = cmd + cufflinkCmd + ' ' + ' '.join(otherParameters) + ' && '
    print cmd[:-3]
    subprocess.call(cmd[:-3],shell=True)
    return result

#===============================================================================
#                 merge cufflinks
#===============================================================================
def merge_cuff_fpkm_cho(path,ConvertFile):
    """This function merges the cufflink fpkm results to one file for CHO genome
    * path: str. Path in which each subfolder is results of one sample.
    * ConvertFile: str. The file that includes 3 columns [geneid, genesymbol,chromosome]
    """
    
    filepath = path #'/data/shangzhong/DE/pgsa/cufflinks'
    os.chdir(filepath)
    folders = [f for f in os.listdir(filepath) if os.path.isdir(os.path.join(filepath,f))]
    folders = natsorted(folders)
    
    # start merge
    data = pd.read_csv(folders[0] +'/genes.fpkm_tracking',sep='\t',header=0,usecols=[0,4,6,9],names=['tracking_id','gene_short_name','locus',folders[0]])
    data['index'] = data['tracking_id'].map(str) + data['gene_short_name'] + data['locus']
    data = data.set_index(['index'])
    for i in range(1,len(folders)):
        df = pd.read_csv(folders[i] +'/genes.fpkm_tracking',sep='\t',header=0,usecols=[0,4,6,9],names=['tracking_id','gene_short_name','locus',folders[i]])
        df['index'] = df['tracking_id'].map(str) + df['gene_short_name'] + df['locus']
        df = df.set_index(['index'])
        df = df[folders[i]]
        data = pd.concat([data,df],axis=1)
     
    data.to_csv(filepath + '/sample_FPKM.csv',sep='\t',index=False)


def merge_cuff_fpkm(path):
    """This function merges cufflinks fpkm results to one file.
    each file has two columns [geneid, genesymbol]
    """
    os.chdir(path)
    folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path,f))]
    folders = natsorted(folders)
    dfs = []
    for f in folders:
        df = pd.read_csv(f+'/genes.fpkm_tracking.txt',sep='\t',header=0,names=['Entrez_ID',f],index_col=0)
        dfs.append(df)
    res_df = pd.concat(dfs,axis=1)
    res_df.to_csv('merged.txt',sep='\t')
#merge_cuff_fpkm('/data/shangzhong/DE/FDA/cufflinks')