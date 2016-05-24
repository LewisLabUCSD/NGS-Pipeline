from ruffus import *
import yaml
from Modules.f01_file_process import *
from Modules.Aligner import BLASR
from Modules.Samtools import *
from Modules.PBhoney import *
import shutil,os
import sys
#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/home/shangzhong/Codes/NewPipeline/Parameters/Pacbio_SV.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
# all parameter
ref_fa = p.ref_fa
# tool parameters
blasr_batch = p.blasr_jobs_per_batch
sam_sort_batch = p.sam_sort_jobs_per_batch

contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('PBHoney start',contact)
os.chdir(file_path)
faFiles = [[os.path.join(file_path,f)] for f in os.listdir(file_path) if f.endswith('.fa')]
Message('Pcbio SV start',contact)
print faFiles

#--------------------- 2. run BLASR -----------------------------------------------------
@jobs_limit(blasr_batch)
@mkdir(faFiles,formatter(),'{path[0]}/sam')
@transform(faFiles,formatter(),'{path[0]}/sam/{basename[0]}.sam')        #regex('.*\.fa'),'.bam')
@check_if_uptodate(check_file_exists)
def run_blasr(input_file,output_file):
    n = num_thread2use(blasr_batch,len(faFiles),thread)
    BLASR(input_file,output_file,ref_fa,n,['-clipping soft'])

#--------------------- 3. Sam2SortBam -----------------------------------------------------
# sam to bam
@follows(run_blasr)
@mkdir(faFiles,formatter(),'{path[0]}/bam')
@transform(run_blasr,formatter('.*\.sam'),'bam/{basename[0]}.bam')
@check_if_uptodate(check_file_exists)
def samtobam(input_file,output_file):
    n = num_thread2use(thread,len(faFiles),thread)
    sam2bam(input_file,output_file,n)
    if os.path.exists('sam'): shutil.rmtree('sam')
# sort bam
@follows(samtobam)
@jobs_limit(sam_sort_batch)
@mkdir(faFiles,formatter(),'{path[0]}/sortBam')
@transform(samtobam,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
@check_if_uptodate(check_file_exists)
def sortbam(input_file,output_file):
    n = num_thread2use(sam_sort_batch,len(faFiles),thread)
    sortBam(input_file,output_file,n)
    if os.path.exists('bam'): shutil.rmtree('bam')
#--------------------- 4. detect SV using PBhoney -----------------------------------------------------
@mkdir(faFiles,formatter(),'{path[0]}/HoneyPie')
@transform(sortbam,formatter('.*\.sort\.bam'),'HoneyPie/{basename[0]}.final.bam')
@check_if_uptodate(check_file_exists)
def Honeypie(input_file,output_file):
    n = num_thread2use(thread,len(faFiles),thread)
    Honey_pie(input_file,output_file,ref_fa,n,'HoneyPie')

@follows(Honeypie)
@mkdir(faFiles,formatter(),'{path[0]}/HoneyTail')
@transform(Honeypie,formatter('.*\.final\.bam'),'HoneyTail/{basename[0]}.tailes')
@check_if_uptodate(check_file_exists)
def Honeytailes(input_file,output_file):
    Honey_tails(input_file,output_file)

@follows(Honeytailes)
@mkdir(faFiles,formatter(),'{path[0]}/HoneySpots')
@transform(Honeypie,formatter('.*\.final\.bam'),'HoneySpots/{basename[0]}.spots')
@check_if_uptodate(check_file_exists)
def Honeyspots(input_file,output_file):
    n = num_thread2use(thread,len(faFiles),thread)
    Honey_spots(input_file,output_file,ref_fa,n)

#---------------------- 5. report succeed -------------------------------------------------------------
@follows(Honeyspots)
def last_function():
    Message('job finished',contact)


if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread, exceptions_terminate_immediately = False,gnu_make_maximal_rebuild_mode = False, 
                 touch_files_only=False,verbose=5)
    except:
        Message('Pacbio SV failed',contact)
        pass
