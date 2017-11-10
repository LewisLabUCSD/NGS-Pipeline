import sarge,sys
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
    