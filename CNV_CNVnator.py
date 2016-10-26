import os,sys
from Modules.f01_file_process import *
from ruffus import *
from Modules.Aligner import bwa_Db,bwa_mem
from Modules.Trimmomatic import Trimmomatic
import yaml
from Modules.Samtools import sortBam
import shutil
from Modules.CNVnator import cnvnator
#============ parameters ======================
#parameter_file =  sys.argv[1]
parameter_file = '/data/shangzhong/Pacbio/CHOS_illu_DNA/cnv/CNVnator.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
# all parameter
ref_fa = p.ref_fa
# trimmomatic parameter
trim = p.trim_reads
trimmomatic = p.trimmomatic_path
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter

bwa_batch = p.bwa_jobs_per_batch
db_path = p.bwa_Db

contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('cnvnator start',contact)
os.chdir(file_path)
fastqFiles = list_fq_files(file_path)
def trim_parameters():
    infiles,outfiles = replace_filename(fastqFiles,'^','trim_')
    for infile, output in zip(infiles,outfiles):
        yield infile,output
#--------------------- 2. trim reads-----------------------------------------------------
@active_if(trim)
@jobs_limit(trim_batch)
@files(trim_parameters)
def trim_reads(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
    remove(input_file)
#--------------------- 3. Map with bwa -----------------------------------------------------
def get_fq():
    fqFiles = list_fq_files(file_path)
    for fq in fqFiles:
        out = 'bam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
        yield fq,out
# build index
@active_if(not os.path.exists(db_path))
@follows(trim_reads)
def bwa_index():
    bwa_Db(db_path,ref_fa)
# align
@jobs_limit(bwa_batch)
@follows(trim_reads,bwa_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/bam')
@files(get_fq)
def run_bwa(input_file,output_file):
    n = num_thread2use(bwa_batch,len(fastqFiles),thread)
    db_index = db_path + '/' + os.listdir(db_path)[0].split('.')[0]
    bwa_mem(input_file,output_file,db_index,n)
#--------------------- 5. Sort bam file --------------------------------------------------
@jobs_limit(trim_batch)
@follows(run_bwa)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
@transform(run_bwa,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
@check_if_uptodate(check_file_exists)
def sort_by_pos(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='pos')
@follows(sort_by_pos)
def remove_bam():
    if os.path.exists('bam'): shutil.rmtree('bam')   # remove bam folder
#------------------- 6. run CNVnator -----------------------------------------------------
@jobs_limit(thread)
@follows(sort_by_pos)
@mkdir(fastqFiles,formatter(),'{path[0]}/cnv')
@transform(sort_by_pos,formatter('.*\.sort\.bam'),'cnv/{basename[0]}.root')
@check_if_uptodate(check_file_exists)
def run_cnvnator(input_file,output_file):
    cnvnator(input_file,output_file)

@follows(run_cnvnator)
def last_function():
    Message('cnvnator succeed',contact)
    
if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
                 touch_files_only=False,verbose=20)
    except:
        Message('cnvnator failed',contact)