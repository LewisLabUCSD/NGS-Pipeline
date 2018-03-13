from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import STAR,STAR_Db
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import sortBam
from Modules.HTseq import htseq_count
from Modules.Homer import *
import yaml
import sys,shutil
from Modules.StringTie import stringtie
#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/data/shangzhong/Proteogenomics/RNAseq_count.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
QC = True
# all parameter
ref_fa = p.ref_fa
annotation = p.gff
gff = p.gff

# fastqc path
#fastqc = p.fastqc_path 
# trimmomatic parameter
trim = p.trim_reads
trim_batch = p.trim_jobs_per_batch
# star parameter
star_batch = p.star_jobs_per_batch
star_db = p.STAR_index_path
other_params = []
run_pass = p.star_pass
contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('fastqc start',contact)
os.chdir(file_path)
fastqFiles = list_fq_files(file_path)
if fastqFiles[0][0].startswith('trim_'):
    trim = False
#--------------------- 2. trim reads-----------------------------------------------------
def trim_parameters():
    fastqFiles = list_fq_files(file_path)
    infiles,outfiles = replace_filename(fastqFiles,'$','.trimmed')
    for infile, output in zip(infiles,outfiles):
        yield infile,output
#------------- run fastqc before trimming -----------
# @active_if(QC)
# @jobs_limit(thread)
    # @mkdir(fastqFiles,formatter(),'{path[0]}/fastqc')
# @files(trim_parameters)
# def run_QC1(input_file,output_file):
#     for fq in input_file:
#         sarge.run('fastqc {input} -o fastqc'.format(input=fq))


#------------ trim file ------------------
@active_if(trim)
#@follows(run_QC1)
@jobs_limit(trim_batch)
@files(trim_parameters)
def trim_reads(input_file,output_file):
    homer_trim_cmd(input_file[0])
    #cmd = 'gzip -c %s.trimmed > trim_%s' % (input_file[0],input_file[0]) #Will have .trimmed appended to file
    #print(cmd)
    #sarge.run(cmd)
    #remove(input_file + '.trimmed')
#------------ run fastqc after trimming ------------
@active_if(QC and trim)
@follows(trim_reads)
@jobs_limit(thread)
@mkdir(fastqFiles,formatter(),'{path[0]}/fastqc')
@files(trim_parameters)
def run_QC2(input_file,output_file):
    for fq in output_file:
    	sarge.run('fastqc {input} -k -o fastqc'.format(input=fq))


# def get_fq():
#     fqFiles = list_fq_files(file_path)
#     for fq in fqFiles:
#         out = 'sortBam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
#         yield fq,out
# # build index
# @active_if(not os.path.exists(star_db))
# @follows(run_QC2)
# def star_index():
#     STAR_Db(star_db,ref_fa,thread)
# # align
# other_params.extend(['--outSAMtype BAM', 'SortedByCoordinate'])
# if run_pass == 2:
#     other_params.append('--twopassMode Basic')
    

# @jobs_limit(star_batch)
# @follows(star_index)
# @mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
# #@transform(fastqFiles,formatter('.*\.f.*?\.gz'),'sortBam/{basename[0]}.bam')
# @files(get_fq)
# def run_star(input_file,output_file):
#     print input_file
#     n = num_thread2use(star_batch,len(fastqFiles),thread)
#     STAR(input_file,output_file,star_db,n,gff,other_params)


@follows(run_QC2)
def last_function():
    Message('fastqc succeed',contact)
    
if __name__ == '__main__':
    try:
        #pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],forcedtorun_tasks = [last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True,
                    touch_files_only=False,verbose=10)#,checksum_level=3)
    except:
        Message('get bam failed',contact)