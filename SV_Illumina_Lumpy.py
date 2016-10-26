import os,sys
from Modules.f01_file_process import *
from ruffus import *
from Modules.Aligner import bwa_Db,bwa_samblaster
from Modules.Trimmomatic import Trimmomatic
import yaml
from Modules.Samtools import sortBam
import shutil
#============ parameters ======================
#parameter_file =  sys.argv[1]
parameter_file = '/home/shangzhong/Codes/NewPipeline/Parameters/SV_lumpy.yaml'
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
bwa_index = p.bwa_index
read_groups = p.read_groups
contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('Lumpy start',contact)
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
@check_if_uptodate(check_file_exists)
def trim_reads(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
    remove(input_file)
#--------------------- 3. Map with bwa -----------------------------------------------------
def get_fq_and_readgroup():
    fqFiles = list_fq_files(file_path)
    for fq, rg in zip(fqFiles,read_groups):
        out = 'bam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
        yield fq,out,rg
# build index
@active_if(not os.path.exists('/'.join(bwa_index.split('/')[:-1])))
@follows(trim_reads)
def bwa_index():
    bwa_Db(bwa_index,ref_fa)
# align
@jobs_limit(bwa_batch)
@follows(bwa_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/bam')
@files(get_fq_and_readgroup)
@check_if_uptodate(check_file_exists)
def run_bwa(input_file,output_file,rg):
    n = num_thread2use(bwa_batch,len(fastqFiles),thread)
    bwa_samblaster(input_file,output_file,bwa_index,n,otherParameters=['-R '+rg+'\\tPL:illumina\\tLB:lib20000\\tPU:unit1'])
#--------------------- 5. Sort bam file --------------------------------------------------
@jobs_limit(trim_batch)
@follows(run_bwa)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
@transform(run_bwa,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
@check_if_uptodate(check_file_exists)
def sort_by_pos(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='pos')
    disc = input_file[:-3] + 'disc.bam'
    disc_sort = disc[:-3] + 'sort.bam'
    sortBam(disc,disc_sort,n,sortType='pos')
    split = input_file[:-3] + 'split.bam'
    split_sort = split[:-3] + 'sort.bam'
    sortBam(input_file,split_sort,n,sortType='pos')
@follows(sort_by_pos)
def remove_bam():
    if os.path.exists('bam'): shutil.rmtree('bam')   # remove bam folder
#--------------------- 6. Get mean and std of paired reads ------------------------------
@jobs_limit(thread)
@follows(sort_by_pos)
@mkdir(fastqFiles,formatter(),'{path[0]/histo}')
@transform(sort_by_pos,formatter('.*\.sort\.bam','histo/{basename[0]}.histo'))
@check_if_uptodate(check_file_exists)
def pairend_histo():






