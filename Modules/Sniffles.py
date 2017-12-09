import sarge,sys

def sniffle(bam,outVCF,otherParameters=['']):
    """run sniffle to detect SV using pacbio"""
    cmd = ('sniffles -m {bam} -v {outVCF} ').format(bam=bam,outVCF=outVCF)
    if otherParameters != ['']:
        cmd = cmd + ' '.join(otherParameters)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

