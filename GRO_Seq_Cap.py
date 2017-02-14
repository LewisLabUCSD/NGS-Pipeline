from ruffus import *
from Modules.f01_file_process import *
from Modules.Trimmomatic import Trimmomatic
from Modules.Aligner import STAR,STAR_Db
from Modules.Samtools import sortBam
from Modules.Homer import *

import yaml
import sys,shutil
#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/data/shangzhong/Proteogenomics/RNAseq_count.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
QC = p.QC
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

contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
# Message('5GRO',contact)
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
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter,22)
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
@active_if(not os.path.exists(db_path))
@follows(trim_reads,run_QC2)
def star_index():
    STAR_Db(db_path,ref_fa,thread)
# align
if trim == False:
    trim_reads=fastqFiles
@jobs_limit(star_batch)
@follows(star_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/bam')
@check_if_uptodate(check_file_exists)
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'bam/{basename[0]}.bam')
def run_star(input_file,output_file):
    n = num_thread2use(star_batch,len(fastqFiles),thread)
    STAR(input_file,output_file,db_path,n,annotation,['--outSAMtype BAM','SortedByCoordinate','--outSAMunmapped Within'])
#--------------------- 4. run Homer ------------------------------------------------------
@follows(run_star)
@mkdir(fastqFiles,formatter(),'{path[0]}/tags')
@check_if_uptodate(check_file_exists)
@transform(run_star,formatter('\.bam'),'tags/{basename[0]}')
def make_tag(input_bam,out_dir):
    make_tag_directory(input_bam,out_dir,ref_fa)

# @follows(run_star)
# def last_function():
#     Message('get bam succeed',contact)
    
if __name__ == '__main__':
    try:
#         pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([run_star],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
                    touch_files_only=False,verbose=15)
    except:
        Message('get bam failed',contact)
    

