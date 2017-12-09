import sarge,sys
def svtyper(in_vcf,out_vcf,bam):
    '''
    this function run svtyper to add genotype to vcf 
    '''
    # 1. generate json file
    json = bam[:-3] + 'json'
    cmd = ('svtyper -B {bam} -l {j} && samtools index {bam}').format(bam=bam,j=json)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    # 2. generate json plot
    sarge.run('lib_stats.R {j} {j}.pdf'.format(j=json))
    # 3. run svtyper
    cmd = ('svtyper -B {bam} -i {invcf} -l {j} -o {out}').format(
                    bam=bam,invcf=in_vcf,j=json,out=out_vcf)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)