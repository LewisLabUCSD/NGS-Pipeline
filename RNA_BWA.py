from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import STAR,STAR_Db, bwa_mem, bwa_Db
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import sortBam
from Modules.HTseq import htseq_count, cufflinks
import yaml
import sys,shutil
#============ parameters ======================
# parameter_file =  sys.argv[1]
import socket
parameter_file =  'Parameters/RNABWA_%s.yaml'%socket.gethostname()
# parameter_file = '/data/shangzhong/Proteogenomics/RNAseq_count.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
QC = p.QC
# all parameter
ref_fa = p.ref_fa.format(use_seq = p.use_seq)
# annotation = p.gff.format(use_seq = p.use_seq)
# trimmomatic parameter
trim = p.trim_reads
trimmomatic = p.trimmomatic_path
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter
# star parameter
star_batch = p.star_jobs_per_batch
db_path = p.STAR_index_path
# htseq parameter
htseq_anno_source = p.htseq_anno_source
strand = p.strand_specific

id_file = p.gene2ref_seq
tax_id = p.tax_id
contact = p.contact
genomeSize = p.genomeSize
sortType = '' # sort by coordinates by default
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('RNA_count start',contact)
os.chdir(file_path)
fastqFiles = list_fq_files(file_path)
if fastqFiles[0][0].startswith('trim_'):
    trim = False
#--------------------- 2. trim reads-----------------------------------------------------
def trim_parameters():
    infiles,outfiles = replace_filename(fastqFiles,'^','trim_')
    for infile, output in zip(infiles,outfiles):
        yield infile,output
#------------- run fastqc before trimming -----------
@active_if(QC)
@jobs_limit(thread)
@mkdir(fastqFiles,formatter(),'{path[0]}/fastqc')
@files(trim_parameters)
def run_QC1(input_file,output_file):
    for fq in input_file:
        sarge.run('fastqc {input} -o fastqc'.format(input=fq))
#------------ trim file ------------------
@active_if(trim)
@follows(run_QC1)
@jobs_limit(trim_batch)
@files(trim_parameters)
def trim_reads(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
    remove(input_file)
#------------ run fastqc after trimming ------------
@active_if(QC and trim)
@follows(trim_reads)
@jobs_limit(thread)
@mkdir(fastqFiles,formatter(),'{path[0]}/fastqc')
@files(trim_parameters)
def run_QC2(input_file,output_file):
    for fq in output_file:
        sarge.run('fastqc {input} -o fastqc'.format(input=fq))
#--------------------- 3. run STAR ------------------------------------------------------
# build index
# @active_if(not os.path.exists(db_path)) ## always active
@follows(trim_reads,run_QC2)
def bwa_index():
    bwa_Db(db_path,ref_fa, genomeSize= genomeSize)
# align
if trim == False:
    trim_reads=fastqFiles
# @jobs_limit(star_batch)
@follows(bwa_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/bam')
# @check_if_uptodate(check_file_exists)
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'bam/{basename[0]}.bam')
def run_bwa(input_file, output_file):

    n = num_thread2use(star_batch,len(fastqFiles),thread)
    bwa_mem(input_file,output_file,
            db_path = os.path.join(db_path,os.path.split(ref_fa)[-1]),
            thread= n)
#--------------------- 4. samtools sort by name -----------------------------------------
@jobs_limit(trim_batch)
@follows(run_bwa)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam_BWA')
# @check_if_uptodate(check_file_exists)
@transform(run_bwa, formatter('.*\.bam'), 'sortBam_BWA/{basename[0]}.sort.bam')
def sort_by_name(input_file,output_file, sortType = sortType):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType=sortType)
    stat = sarge.get_stdout('samtools flagstat {bam}'.format(bam=output_file))
    with open(output_file[:-3]+'flagstat.txt','w') as f:
        f.write(stat)
@follows(sort_by_name)
def remove_bam():
    # pass
    if os.path.exists('bam'): shutil.rmtree('bam')   # remove bam folder
#--------------------- 5.2 run cufflinks---------------------------------------------------

@active_if(p.cufflinks)
@jobs_limit(thread)
@follows(remove_bam)
@mkdir(fastqFiles,formatter(),'{path[0]}/cufflinks_BWA')
# @check_if_uptodate(check_file_exists)
@transform(sort_by_name,formatter('.*\.sort\.bam'),'cufflinks_BWA/{basename[0]}')
def run_cufflinks(input_file,output_file):
    cufflinks(input_file,output_file,annotation,strand)


#--------------------- 5.2 run htseq -----------------------------------------------------
@active_if(not p.cufflinks)
@follows(remove_bam)
@mkdir(fastqFiles,formatter(),'{path[0]}/htseq_BWA')
@check_if_uptodate(check_file_exists)
@transform(sort_by_name,formatter('.*\.sort\.bam'),'htseq_BWA/{basename[0]}.txt')
def run_htseq(input_file,output_file):
    htseq_count(input_file,output_file,annotation,strand,htseq_anno_source)
#--------------------- 6. ID convertion -----------------------------------------------------
@active_if((htseq_anno_source=='ncbi') & (not p.cufflinks))
@follows(run_htseq)
@transform(run_htseq,suffix('.txt'),'.count.txt')
def id_convert(input_file,output_file):
    print(input_file+ '--->' + output_file)
    id_symbol_conversion(input_file,output_file,id_file,tax_id)
#--------------------- 7. return finish message -----------------------------------------------------
if (htseq_anno_source == 'ncbi') & (not p.cufflinks) :
    @follows(run_htseq,id_convert)
    def last_function():
        Message('RNA_count finished',contact)
elif p.cufflinks:
    @follows(run_cufflinks)
    def last_function():
        Message('RNA_count finished',contact)
else:
    @follows(run_htseq)
    def last_function():
        Message('RNA_count finished',contact)

if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
                 touch_files_only=False,verbose=5)
    except:
        Message('RNA_count failed',contact)
    