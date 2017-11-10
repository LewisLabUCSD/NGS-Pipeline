import sarge,sys,os
def lumpyexpress(in_bams,out_vcf,others=['']):
    '''This function runs lumpy express
    * in_bams: sorted bams '''
    bams = ','.join(in_bams)
    splits = ','.join([b[:-3]+'split.bam' for b in in_bams])
    discs = ','.join([b[:-3]+'disc.bam' for b in in_bams])
    cmd = ('lumpyexpress -B {bams} -S {splits} -D {discs} {other} '
              '-k -o {out}').format(
                bams=bams,splits=splits,discs=discs,out=out_vcf,other=' '.join(others))
    with open('cmd.sh','w') as f:
        f.write('#!/bin/bash\n' + cmd)
    print(cmd);sys.stdout.flush()
    sarge.run('chmod 777 cmd.sh && ./cmd.sh && rm cmd.sh')
    