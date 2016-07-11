from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import STAR,STAR_Db
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import sortBam
from Modules.HTseq import htseq_count
import yaml
import sys,shutil

#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/data/shangzhong/DE/winzeler/RNAseq_count.yaml'
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
adapter = p.adapter
# star parameter
star_batch = p.star_jobs_per_batch
db_path = p.STAR_index_path
# htseq parameter
htseq_path = p.htseq_path
htseq_anno_source = p.htseq_anno_source
strand = p.strand_specific

id_file = p.gene2ref_seq
tax_id = p.tax_id
contact = p.contact
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

@active_if(trim)
@jobs_limit(trim_batch)
@files(trim_parameters)
def trim_reads(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
    remove(input_file)
#--------------------- 3. run STAR ------------------------------------------------------
# build index
@active_if(not os.path.exists(db_path))
@follows(trim_reads)
def star_index():
    STAR_Db(db_path,ref_fa,thread)
# align
if trim == True:
    @jobs_limit(star_batch)
    @follows(trim_reads,star_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/bam')
    @check_if_uptodate(check_file_exists)
    @transform(trim_reads,formatter('.*\.f.*?\.gz'),'bam/{basename[0]}.bam')
    def run_star(input_file,output_file):
        n = num_thread2use(star_batch,len(fastqFiles),thread)
        STAR(input_file,output_file,db_path,n,annotation,['--outSAMtype BAM','Unsorted'])
else:
    @jobs_limit(star_batch)
    @follows(star_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/bam')
    @check_if_uptodate(check_file_exists)
    @transform(fastqFiles,formatter('.*\.f.*?\.gz'),'bam/{basename[0]}.bam')
    def run_star(input_file,output_file):
        n = num_thread2use(star_batch,len(fastqFiles),thread)
        STAR(input_file,output_file,db_path,n,annotation,['--outSAMtype BAM','Unsorted'])
#--------------------- 4. samtools sort by name -----------------------------------------
@jobs_limit(trim_batch)
@follows(run_star)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
@check_if_uptodate(check_file_exists)
@transform(run_star,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
def sort_by_name(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='name')
@follows(sort_by_name)
def remove_bam():
    if os.path.exists('bam'): shutil.rmtree('bam')   # remove bam folder
#--------------------- 5. run htseq -----------------------------------------------------
@follows(remove_bam)
@mkdir(fastqFiles,formatter(),'{path[0]}/htseq')
@check_if_uptodate(check_file_exists)
@transform(sort_by_name,formatter('.*\.sort\.bam'),'htseq/{basename[0]}.txt')
def run_htseq(input_file,output_file):
    htseq_count(input_file,output_file,annotation,strand,htseq_path,htseq_anno_source)
#--------------------- 6. ID convertion -----------------------------------------------------
@active_if(htseq_anno_source=='ncbi')
@follows(run_htseq)
@transform(run_htseq,suffix('.txt'),'.count.txt')
def id_convert(input_file,output_file):
    print(input_file+ '--->' + output_file)
    id_symbol_conversion(input_file,output_file,id_file,tax_id)
#--------------------- 7. return finish message -----------------------------------------------------
if htseq_anno_source == 'ncbi':
    @follows(run_htseq,id_convert)
    def last_function():
        Message('RNA_count finished',contact)
else:
    @follows(run_htseq)
    def last_function():
        Message('RNA_count finished',contact)

if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = False, 
                 touch_files_only=False)
    except:
        Message('RNA_count failed',contact)
    