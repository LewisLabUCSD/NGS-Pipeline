import os
import sarge

def htseq_count(sortedBam,countFile,annotation,strand,outpath,annotationSource):
    """This function run htseq_count to count reads given bam file
    * sortedBam: str. Bamfile name
    * countFile: outputfilename
    * annotation: annotation file
    * outputpath: path to store the result files
    * annotation: source. 'ncbi','ensembl'
    """
    if not os.path.exists(outpath):(outpath)
    # 1. check whether outputpath exist
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    # 2. check the annotation source
    if annotationSource == 'ncbi':
        seqType = 'exon'
        id_attr = 'gene'
    elif annotationSource == 'ensembl':
        seqType = 'exon'
        id_attr = 'gene_id'
    elif annotationSource == 'genedb':
        seqType = 'CDS'
        id_attr = 'Parent'
    # 3. run htseq-count
    cmd = ('htseq-count -f bam -s {strand} -t {type} -i {gene} {bam} {annotation} > {output}').format(strand=strand,
         type=seqType,gene=id_attr,bam=sortedBam,annotation=annotation,output=countFile)#os.path.join(outpath,countFile))
    print(cmd)
    sarge.run(cmd)
        
        
def Message(string,email):
    """
    This function send message to email when it run. 
    Used to calculate the time code runs.
    """
    cmd = ('echo {quote}|mailx -s "{string}" {email}').format(quote="",string=string,email=email)
    sarge.run(cmd)
    
