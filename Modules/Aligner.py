import sarge
import os
import shutil
import sys
#===============================================================================
#                         STAR
#===============================================================================
def STAR_Db(db_path,ref_fa,thread=1,annotation = '',genomeSize='large'):
    """
    This function generates database for alignment using STAR
    """
    if not os.path.exists(db_path): os.mkdir(db_path)
    if os.listdir(db_path) == []:
        cmd = ('STAR --runMode genomeGenerate --genomeDir {db_path} '
               '--genomeFastaFiles {ref_fa} --runThreadN {thread} '
               '--limitGenomeGenerateRAM 100000000000 ').format(
                db_path=db_path,ref_fa=ref_fa,thread=str(thread))
        if annotation != '':
            cmd = cmd + ('--sjdbGTFfile {gff3} --sjdbGTFtagExonParentTranscript Parent '
                         '--sjdbOverhang 100').format(gff3=annotation)   # for geneDb add --sjdbGTFfeatureExon CDS
        if genomeSize == 'small':
            cmd = cmd + '--genomeChrBinNbits 6 --genomeSAindexNbases 4'
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

def STAR(fastqFiles,outSamFile,db_path,thread=1,annotation='',otherParameters=['']):
    """STAR for single end read"""
    if annotation != '':
        otherParameters.extend(['--sjdbGTFfile {gff}'.format(gff=annotation)])
    if annotation.endswith('gff') or annotation.endswith('gff3'):
        otherParameters.append('--sjdbGTFtagExonParentTranscript Parent')
    # generate command
    if len(fastqFiles) == 1:
        starCmd = ('STAR --genomeDir {ref} --readFilesCommand zcat '
                     '--readFilesIn {fq1} --runThreadN {thread} '
                     '--outFileNamePrefix {output} --outSAMstrandField intronMotif '
                     '--outFilterIntronMotifs RemoveNoncanonical').format(
                    ref=db_path,fq1=fastqFiles[0],
                    thread=thread,output=outSamFile)
    elif len(fastqFiles) == 2:
        starCmd = ('STAR --genomeDir {ref} --readFilesCommand zcat '
                     '--readFilesIn {fq1} {fq2} --runThreadN {thread} '
                     '--outFileNamePrefix {output} --outSAMstrandField intronMotif '
                     '--outFilterIntronMotifs RemoveNoncanonical').format(
                    ref=db_path,fq1=fastqFiles[0],fq2=fastqFiles[1],
                    thread=thread,output=outSamFile)
    cmd = starCmd + ' ' + ' '.join(otherParameters)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    if 'SortedByCoordinate' in otherParameters:
        outFile = outSamFile+'Aligned.sortedByCoord.out.bam'
    else:
        outFile = outSamFile+'Aligned.out.bam'
    os.rename(outFile,outSamFile)
    if os.path.exists(outSamFile+'_STARgenome'):
        shutil.rmtree(outSamFile+'_STARgenome')
    
    
def BLASR(faFile,outBam,ref_fa,thread,otherParameters=['']):
    """This function runs BLASR"""
    
    if otherParameters != ['']:
        other = ' '.join(otherParameters)
    else:
        other = ''
    cmd = ('blasr {input} {ref} -sam -nproc {thread} {other} | samtools view -hb - > {out}').format(
                    input=faFile,ref=ref_fa,thread=str(thread),other=other,out=outBam,)
    
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

#===============================================================================
#                         bwa
#===============================================================================
def bwa_Db(db_path,ref_fa):
    """build bwa index"""
    if not os.path.exists(db_path):
        os.mkdir(db_path)
    cmd = ('bwa index -p {db_path}/bwa -a bwtsw {fa}').format(fa=ref_fa,db_path=db_path)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    
    
def bwa_mem(fqFile,outSam,db_name,thread,otherParameters=['']):
    """run bwa"""
    if otherParameters != ['']:
        other =  ' '.join(otherParameters) + ' '
    else:
        other = ''
    if len(fqFile) == 1:
        bwaCmd = ('bwa mem -t {thread} {other}{db} {fq} | samtools view -bh - > {out} ').format(
                    thread=str(thread),other=other,db=db_name,fq=fqFile[0],
                    out=outSam)
    else:
        bwaCmd = ('bwa mem -t {thread} {other}{db} {fq1} {fq2} | samtools view -bh - > '
        '{out} ').format(thread=str(thread),other=other,db=db_name,fq1=fqFile[0],
        fq2=fqFile[1],out=outSam)
    print(bwaCmd);sys.stdout.flush()
    sarge.run(bwaCmd)
#bwa_mem('/data/shangzhong/Pacbio/sniffle/CHOS.fq.gz','/data/shangzhong/Pacbio/sniffle/result.bam','5',['-x pacbio'])   


def bwa_samblaster(fqFiles,outBam,db_name,thread,otherParameters=['']):
    '''map for lumpy '''
    if len(fqFiles) != 2:
        assert False,'fastq files are not paired'
    if otherParameters != ['']:
        other =  ' '.join(otherParameters) + ' '
    else:
        other = ''
    split = outBam[:-3]+'split.sam'
    disc = outBam[:-3] + 'disc.sam'
    cmd = ('bwa mem -t {thread} {other}{db} {fq1} {fq2} | samblaster --addMateTags -e -d {disc} -s {split} | \
            samtools view -Sb - > {out}').format(thread=str(thread),other=other,db=db_name,fq1=fqFiles[0],
            fq2=fqFiles[1],disc=disc,split=split,out=outBam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

#===============================================================================
#                         HISAT2
#===============================================================================
def hisat2_Db(ref_fa,db,thread=1):
    """
    """
    cmd = ('hisat2-build -p {t} {ref} {name} ').format(t=str(thread),ref=ref_fa,name=db)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    


def hisat2(fqFile,outBam,db_name,thread,otherParameters=['']):
    """
    """
    if otherParameters != ['']:
        other =  ' '.join(otherParameters) + ' '
    else:
        other = ''
    if len(fqFile) == 1:
        hisat2Cmd = ('hisat2 -x {db} -U {fq} -t {other} -p {thread} '
                     '| samtools view -bh - > {out}').format(db=db_name,fq=fqFile[0],
                      other=other,thread=str(thread),out=outBam)
    else:
        hisat2Cmd = ('hisat2 -x {db} -1 {fq1} -2 {fq2} -t {other} -p {thread} '
                     '| samtools view -bh - > {out}').format(db=db_name,fq1=fqFile[0],fq2=fqFile[1],
                    other=other,thread=str(thread),out=outBam)
    
    print(hisat2Cmd);sys.stdout.flush()
    sarge.run(hisat2Cmd)

#===============================================================================
#                     ngmlr
#===============================================================================
def ngmlr(in_fa,outBam,ref_fa,thread):
    '''run nglmr for better SV detection using pacbio'''
    cmd = ('ngmlr -t {thread} -r {ref} -q {fa} | samtools view -hb - > outBam').format(
                            thread=str(thread),ref=ref_fa,fa=in_fa)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)