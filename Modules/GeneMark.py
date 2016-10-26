import sarge,os,sys
def geneMark_ES(ref_fa,other_params=['']):
    '''run geneMark_ES'''
    cmd = ('gmes_petap.pl --ES {other} --sequence {fa}').format(fa=ref_fa,
                                            other=' '.join(other_params))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    #return os.getcwd() +'/genemark.gff'