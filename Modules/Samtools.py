import sarge

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
    print(cmd)
    sarge.run(cmd)
    if sortType !='name':
        cmd = ('samtools index {bam} ').format(bam=sortedBamFile)
        print(cmd)
        sarge.run(cmd)

def sam2bam(samFile,bamFile,thread):
    """
    This function change sam file to bam file
    """
    cmd = ('samtools view -@ {thread} -h {sam} -o {bam} ').format(
            thread=thread,sam=samFile,bam=bamFile)
    print(cmd)
    sarge.run(cmd)
            
