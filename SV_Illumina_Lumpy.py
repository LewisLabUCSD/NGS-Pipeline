import os,sys
from Modules.f01_file_process import *
from ruffus import *
from Modules.Aligner import bwa_Db,bwa_samblaster
from Modules.Trimmomatic import Trimmomatic
import yaml
from Modules.Samtools import sortBam,merge_bams
import shutil
from Modules.Lumpy import lumpyexpress
from Modules.SVTyper import svtyper
#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/data/shangzhong/SV_lumpy/SV_lumpy.yaml'
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
def run_bwa_index():
    bwa_Db(bwa_index,ref_fa)
# align
@jobs_limit(bwa_batch)
@follows(trim_reads,run_bwa_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/bam')
@files(get_fq_and_readgroup)
def run_bwa(input_file,output_file,rg):
    n = num_thread2use(bwa_batch,len(fastqFiles),thread)
    lib = rg.split('\\t')[2][3:]
    readgroup = '\'' + rg+'\\tLB:'+lib+'\\tPL:illumina\\tPU:unit1\''
    bwa_samblaster(input_file,output_file,bwa_index,n,otherParameters=['-R '+ readgroup])
#--------------------- 5. Sort bam file --------------------------------------------------
@jobs_limit(trim_batch)
@follows(run_bwa)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
@transform(run_bwa,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
@check_if_uptodate(check_file_exists)
def sort_by_pos(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='pos')
    disc = input_file[:-3] + 'disc.sam'
    disc_sort = output_file[:-3] + 'disc.bam'
    sortBam(disc,disc_sort,n,sortType='pos')
    split = input_file[:-3] + 'split.sam'
    split_sort = output_file[:-3] + 'split.bam'
    sortBam(split,split_sort,n,sortType='pos')
@follows(sort_by_pos)
def remove_bam():
    if os.path.exists('bam'): shutil.rmtree('bam')   # remove bam folder
#--------------------- 6. run lumpyexpress ------------------------------
@follows(sort_by_pos)
@mkdir(fastqFiles,formatter(),'{path[0]}/vcf')
@merge(sort_by_pos,'vcf/lumpy.vcf')
@check_if_uptodate(check_file_exists)
def run_lumpyexpress(input_files,output_file):
    lumpyexpress(input_files,output_file,['-T vcf'])
#--------------------- 7. run SVTyper ------------------------------
@follows(run_lumpyexpress)
@mkdir(fastqFiles,formatter(),'{path[0]}/merge')
@merge(sort_by_pos,'merge/merge.bam')
@check_if_uptodate(check_file_exists)
def run_svtyper(input_files,output_file):
    if len(input_files) > 1:
        merge_bams(input_files,output_file)
    else:
        os.rename(input_files[0],output_file)
    sarge.run('samtools index {b}'.format(b=output_file))
    svtyper('vcf/lumpy.vcf','vcf/lumpy_gt.vcf',output_file)
    shutil.move('merge/merge.json','vcf')
    shutil.move('merge/merge.json.pdf','vcf')

@follows(run_svtyper)
def last_function():
    Message('lumpy succeed',contact)
    
if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = False, 
                 touch_files_only=False,verbose=20)
    except:
        Message('lumpyexpress failed',contact)




