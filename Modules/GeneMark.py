import pandas as pd
import sarge,os,sys
def geneMark_ES(ref_fa,other_params=['']):
    '''run geneMark_ES'''
    cmd = ('gmes_petap.pl --ES {other} --sequence {fa}').format(fa=ref_fa,
                                            other=' '.join(other_params))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    sarge.run('genemark_gtf2gff3 genemark.gtf > genemark.gff')  # this code is download from maker
    df = pd.read_csv('genemark.gff',sep='\t',comment='#',header=None)
    df[0] = df[0].map(lambda x: x.split(' ')[0])
    df.to_csv('genemark.gff3',sep='\t',index=False,header=None)
    #return os.getcwd() +'/genemark.gff'

