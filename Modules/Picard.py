import sarge
import sys,os

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
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

# if __name__ == '__main__':
#     path = '/data/shangzhong/DetectVirus/unmap_bam'
#     picard = '/home/shangzhong/Installation/picard-tools-1.141/picard.jar'
#     os.chdir(path)
#     bams = [f for f in os.listdir(path) if f.endswith('.bam')]
#     for bam in bams:
#         out = bam.split('.')[0]
#         sam2fq(bam,out,picard,'single')


def build_fa_dict(ref_fa,picard):
    '''build dictionary file for fa file '''
    out = '.'.join(ref_fa.split('.')[:-1]) + '.dict'
    cmd = ('java -jar {picard} CreateSequenceDictionary R={ref} OO={out}').format(
            picard = picard,ref=ref_fa,out=out)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def mark_duplicates(sortBam,dedupBam,picard):
    '''mark duplicates'''
    cmd = ('java -Djava.io.tmpdir=tmp -jar {picard} MarkDuplicates I={input} O={out} '
           'CREATE_INDEX=true METRICS_FILE=metrics.txt MAX_RECORDS_IN_RAM=8000000 '
           'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 '
           'VALIDATION_STRINGENCY=LENIENT').format(picard=picard,input=sortBam,out=dedupBam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def add_readgroup(sortBam,rgBam,readgroup,picard):
    '''add read group'''
    if not os.path.exists('tmp'):os.mkdir('tmp')
    rg = readgroup.split('\\t')
    ID = rg[1][3:]
    SM = rg[2][3:]
    PL = 'illumina'
    LB = 'lib20000'
    PU = 'unit1'
    cmd = ('java -jar {picard} AddOrReplaceReadGroups I={input} O={rgBam} SO=coordinate '
            'RGID={ID} RGSM={SM} RGPL={PL} RGLB={LB} RGPU={PU} TMP_DIR=tmp').format(
           picard=picard,input=sortBam,rgBam=rgBam,ID=ID,SM=SM,PL=PL,LB=LB,PU=PU)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

    
    