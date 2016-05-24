import sarge

def sam2fq(inSam,outPrefix,picard,endType):
    """sam file to fastq files
    """
    if endType == 'single':
        cmd = ('java -jar {picard} SamToFastq I={input} F={fq} '
               'VALIDATION_STRINGENCY=LENIENT ').format(
            picard=picard,input=inSam,fq=outPrefix+'.fq.gz')
    else:
        cmd = ('java -jar {picard} SamToFastq I={input} F={fq1} F2={fq2} '
               'VALIDATION_STRINGENCY=LENIENT ').format(picard=picard,
                input=inSam,fq1=outPrefix+'_1.fq.gz',fq2=outPrefix+'_2.fq.gz')
    sarge.run(cmd)

# if __name__ == '__main__':
#     path = '/data/shangzhong/DetectVirus/unmap_bam'
#     picard = '/home/shangzhong/Installation/picard-tools-1.141/picard.jar'
#     os.chdir(path)
#     bams = [f for f in os.listdir(path) if f.endswith('.bam')]
#     for bam in bams:
#         out = bam.split('.')[0]
#         sam2fq(bam,out,picard,'single')