import os,sys
from Modules.f01_file_process import *
from ruffus import *
from Modules.Aligner import bwa_mem,bwa_Db,ngmlr
from Modules.Sniffles import sniffle
import yaml
from Modules.Samtools import sortBam

#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/home/shangzhong/Codes/NewPipeline/Parameters/SV_Pacbio_Sniffle.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
# all parameter
ref_fa = p.ref_fa
db_path = p.bwa_db
contact = p.contact
aligner = p.aligner
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('Sniffle start',contact)
os.chdir(file_path)
#--------------------- 2. align all files -----------------------------------------------
if aligner == 'bwa':
    fastqFiles = [f for f in os.listdir(file_path) if f.endswith('fq.gz') or f.endswith('fastq.gz')]
    # build index
    @active_if(not os.path.exists(db_path))
    def bwa_index():
        bwa_Db(db_path,ref_fa)
        os.chdir(file_path)
     
    @follows(bwa_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/bam')
    @check_if_uptodate(check_file_exists)
    @transform(fastqFiles,formatter('.*\.f.*q\.gz'),'bam/{basename[0]}.bam')
    def run_bwa(input_file,output_file):
        print(input_file + '-->' + output_file)
        bwa_mem([input_file],output_file,db_path+'/bwa',thread,otherParameters=['-M','-x pacbio'])
elif aligner == 'ngmlr':
    faFiles = [f for f in os.listdir(file_path) if f.endswith('fa.gz') or f.endswith('fasta.gz')]
    @transform(faFiles,formatter('.*\.f.*a\.gz'),'bam/{basename[0]}.bam')
    def run_ngmlr(input_file,output_file):
        ngmlr(input_file,output_file,ref_fa,thread)
#--------------------- 3. sort bam file -----------------------------------------------
@follows(run_bwa,run_ngmlr)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
@check_if_uptodate(check_file_exists)
@transform(run_bwa,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
def run_sortBam(input_file,output_file):
    sortBam(input_file,output_file,thread)
#--------------------- 4. Detect SV -----------------------------------------------
@follows(run_sortBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/vcf')
@check_if_uptodate(check_file_exists)
@transform(run_sortBam,formatter('.*\.sort\.bam'),'vcf/{basename[0]}.vcf')
def run_sniffle(input_file,output_file):
    sniffle(input_file,output_file,otherParameters=[''])
#--------------------- 4. return finish message -----------------------------------------------------
@follows(run_sniffle)
def last_function():
    Message('SV_Sniffle finished',contact)
    

if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
                 touch_files_only=False,verbose=15)
    except:
        Message('SV_Sniffle failed',contact)