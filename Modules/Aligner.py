import sarge
import os
import shutil

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
    print(cmd)
    sarge.run(cmd)

def STAR(fastqFiles,outSamFile,db_path,thread=1,annotation='',otherParameters=['']):
    """STAR for single end read"""
    if annotation != '':
        otherParameters.extend(['--sjdbGTFfile {gff}'.format(gff=annotation)])
    if annotation.endswith('gff'):
        otherParameters.append('--sjdbGTFtagExonParentTranscript Parent')
    # generate command
    if len(fastqFiles) == 1:
        starCmd = ('STAR --genomeDir {ref} --readFilesCommand zcat '
                     '--readFilesIn {fq1} --runThreadN {thread} '
                     '--outFileNamePrefix {output} ').format(
                    ref=db_path,fq=fastqFiles[0],
                    thread=thread,output=outSamFile)
    elif len(fastqFiles) == 2:
        starCmd = ('STAR --genomeDir {ref} --readFilesCommand zcat '
                     '--readFilesIn {fq1} {fq2} --runThreadN {thread} '
                     '--outFileNamePrefix {output} '
                     ).format(
                    ref=db_path,fq1=fastqFiles[0],fq2=fastqFiles[1],
                    thread=thread,output=outSamFile)
    cmd = starCmd + ' ' + ' '.join(otherParameters)
    print(cmd)
    sarge.run(cmd)
    os.rename(outSamFile+'Aligned.out.bam',outSamFile)
    shutil.rmtree(outSamFile+'_STARgenome')
    
    
def BLASR(faFile,outSam,ref_fa,thread,otherParameters=['']):
    """This function runs BLASR"""
    
    cmd = ('blasr {input} {ref} -sam -out {out} -nproc {thread} ').format(
                    input=faFile,ref=ref_fa,out=outSam,thread=str(thread))
    if otherParameters != ['']:
        cmd = cmd + ' '.join(otherParameters)
    print(cmd)
    sarge.run(cmd)
    