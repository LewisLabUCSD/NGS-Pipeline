import sarge,sys,glob
import pandas as pd

def stringtie(in_bam,out_gtf,thread,annotation):
    '''
    '''
    quant = out_gtf[:-3] + 'abund.tab'
    cov_ref = out_gtf[:-3] + 'cov_ref.gtf'
    cmd = ('stringtie {bam} -o {gtf} -p {t} -G {gff} -A {q} \
            -C {cov}').format(bam=in_bam,gtf=out_gtf,t=str(thread),
                                    gff=annotation,q=quant,cov=cov_ref)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def merge_stringtie_tpm(path):
    """This function merges tpm fpkm results to one file.
    each file has two columns [geneid, tpm]
    """
    files = glob.glob(path + '/*.tab')
    dfs = []
    for f in files:
        sp = f.split('/')[-1].split('.')[0]
        df = pd.read_csv(f,sep='\t',header=0,usecols=[1,8],names=['name',sp],index_col=0)
        df = df[df.index.values !='-']
        df = df.groupby('name').sum()
        dfs.append(df)
    res_df = pd.concat(dfs,axis=1)
    return res_df