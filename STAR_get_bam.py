from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import STAR_Db,STAR
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import build_fa_index,merge_bams
from Modules.Picard import build_fa_dict,mark_duplicates,add_readgroup
from Modules.GATK import *
import yaml
import shutil
import glob


#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/home/shangzhong/Codes/NewPipeline/Parameters/STAR_get_bam.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
# all parameter
ref_fa = p.ref_fa
gff = p.gff
# trimmomatic parameter
trim = p.trim_reads
trimmomatic = p.trimmomatic_path
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter

picard = p.picard

star_batch = p.star_jobs_per_batch
star_db = p.star_index


contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#Message('get bam start',contact)
os.chdir(file_path)
#===============================================================================
#                     Part I. Preprocess
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
fastqFiles = list_fq_files(file_path)
if fastqFiles[0][0].startswith('trim_'):
    trim = False
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
#--------------------- 4. Map with STAR -----------------------------------------------------
def get_fq():
    fqFiles = list_fq_files(file_path)
    for fq in fqFiles:
        out = 'sortBam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
        yield fq,out
# build index
@active_if(not os.path.exists(star_db))
@follows(trim_reads)
def star_index():
    STAR_Db(star_db,ref_fa,thread)
# align
if gff == '':
    @jobs_limit(star_batch)
    @follows(trim_reads,star_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
    @check_if_uptodate(check_file_exists)
    @files(get_fq)
#     @transform(trim_reads,formatter('.*\.f.*?\.gz'),'f01_bam/{basename[0]}.sort.bam')
    def run_star(input_file,output_file):
        n = num_thread2use(star_batch,len(fastqFiles),thread)
        STAR(input_file,output_file,star_db,n,'',['--outSAMtype BAM','SortedByCoordinate','--twopassMode Basic'])
else:
    @jobs_limit(star_batch)
    @follows(trim_reads,star_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
    @check_if_uptodate(check_file_exists)
    @files(get_fq)
#     @transform(fastqFiles,formatter('.*\.f.*?\.gz'),'f01_bam/{basename[0]}.sort.bam')
    def run_star(input_file,output_file):
        n = num_thread2use(star_batch,len(fastqFiles),thread)
        STAR(input_file,output_file,star_db,n,'',['--outSAMtype BAM','SortedByCoordinate'])
        
@follows(run_star)
def last_function():
    Message('get bam succeed',contact)
    
if __name__ == '__main__':
    try:
#         pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
        touch_files_only=False)
    except:
        Message('get bam failed',contact)