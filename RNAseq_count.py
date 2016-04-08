from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import STAR,STAR_Db
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import sortBam
from Modules.HTseq import htseq_count
import yaml
import sys
#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/home/shangzhong/Codes/NewPipeline/Parameters/RNAseq_count.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
# all parameter
ref_fa = p.ref_fa
annotation = p.gff
# trimmomatic parameter
trim = p.trim_reads
trimmomatic = p.trimmomatic_path
trim_batch = p.trim_jobs_per_batch
phred = p.phred
adapter = p.adapter
# star parameter
star_batch = p.star_jobs_per_batch
db_path = p.STAR_index_path
# htseq parameter
htseq_path = p.htseq_path
htseq_anno_source = p.htseq_anno_source
strand = p.strand_specific

id_file = p.id_symbol_convert_file
contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
os.chdir(file_path)
fastqFiles = list_fq_files(file_path)
#--------------------- 2. trim reads-----------------------------------------------------
def trim_parameters():
    infiles,outfiles = replace_filename(fastqFiles,'^','trim_')
    for infile, output in zip(infiles,outfiles):
        yield infile,output

@active_if(trim)
@jobs_limit(trim_batch)
@files(trim_parameters)
def trim_reads(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    Trimmomatic(input_file,output_file,trimmomatic,n,phred,adapter)
    remove(input_file)
#--------------------- 3. run STAR ------------------------------------------------------
# build index
@active_if(not os.path.exists(db_path))
@follows(trim_reads)
def star_index():
    STAR_Db(db_path,ref_fa,thread)
# align
n = num_thread2use(trim_batch,len(fastqFiles),thread)
if trim:
    @follows(trim_reads,star_index)
    @transform(trim_reads,regex('\.f.*?\.gz'),'.bam')
    @jobs_limit(star_batch)
    def run_star(input_file,output_file):
        STAR(input_file,output_file,db_path,n,annotation,['--outSAMtype BAM','Unsorted'])
else:
    @follows(trim_reads,star_index)
    @transform(fastqFiles,regex('\.f.*?\.gz'),'.bam')
    @jobs_limit(star_batch)
    def run_star(input_file,output_file):
        STAR(input_file,output_file,db_path,n,annotation,['--outSAMtype BAM','Unsorted'])
#--------------------- 4. samtools sort by name -----------------------------------------
@jobs_limit(trim_batch)
@transform(run_star,suffix('.bam'),'.sort.bam')
def sort_by_name(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='name')
    remove(input_file)
#--------------------- 5. run htseq -----------------------------------------------------
@follows(sort_by_name,mkdir(htseq_path))
@transform(sort_by_name,formatter(".+sort.bam"),htseq_path+"/{basename[0]}.txt")
def run_htseq(input_file,output_file):
    htseq_count(input_file,output_file,annotation,strand,htseq_path,htseq_anno_source)
    remove(input_file)
#--------------------- 6. ID convertion -----------------------------------------------------
@active_if(htseq_anno_source=='ncbi')
@transform(run_htseq,suffix('.txt'),'.count.txt')
def id_convert(input_file,output_file):
    print input_file+ '--->' + output_file
    id_symbol_conversion(input_file,output_file,id_file)
#--------------------- 7. return finish message -----------------------------------------------------
if htseq_anno_source == 'ncbi':
    @follows(id_convert)
    def last_function():
        Message('job finished',contact)
else:
    @follows(run_htseq)
    def last_function():
        Message('job finished',contact)

if __name__ == '__main__':
    try:
        pipeline_run(multiprocess=thread, exceptions_terminate_immediately = True,gnu_make_maximal_rebuild_mode = False, 
                 touch_files_only=False)
    except:
        Message('job failed',contact)
    