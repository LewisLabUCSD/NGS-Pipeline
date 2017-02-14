import sarge
import sys

def sortBam(bamFile,sortedBamFile,thread=1,sortType=''):
    """
    This function sort bam files
    """
    if sortType == 'name':
            tag = ' -n'
    else:
        tag = ''
    cmd = ('samtools sort{tag} -m 4G -@ {thread} -T {sort} -o {sortBam} {bam} ').format(
            tag=tag,thread=str(thread),sort=bamFile[:-3]+'sort',bam=bamFile,sortBam=sortedBamFile)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    if sortType !='name':
        cmd = ('samtools index {bam} ').format(bam=sortedBamFile)
        print(cmd);sys.stdout.flush()
        sarge.run(cmd)

def sam2bam(samFile,bamFile,thread):
    """
    This function change sam file to bam file
    """
    cmd = ('samtools view -@ {thread} -h {sam} -o {bam} ').format(
            thread=thread,sam=samFile,bam=bamFile)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
            

def build_fa_index(ref_fa):
    '''build fai file for fa file for GATK
    '''
    cmd = ('samtools faidx {ref}').format(ref=ref_fa)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def merge_bams(bamfiles,outputbam):
    """this function merges bam files into one"""
    if len(bamfiles) == 1:
        cmd = ('mv {input} {output}').format(input=bamfiles[0],output=outputbam)
    else:
        bam = ' '.join(bamfiles)
        cmd = ('samtools merge -f {output} {input}').format(output=outputbam,input=bam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
