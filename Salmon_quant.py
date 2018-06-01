from ruffus import *
from Modules.f01_file_process import *
from Modules.Trimmomatic import conda_Trimmomatic
from Modules.Salmon import *

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
gff_fn = p.gff
gff_type = p.gff_type
salmonDb = p.salmon_index
# trimmomatic parameter
trim = p.trim_reads
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter
# star parameter
lib_type = p.lib_type
salmon_batch = p.salmon_batch
contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
# Message('salmon start',contact)
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
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
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
#------------------- 3. run salmon ------------------
@active_if(not os.path.exists(salmonDb))
@follows(run_QC2,trim_reads)
@jobs_limit(thread)
def run_salmon_index():
    salmon_index(salmonDb,ref_fa)

if trim == False:
    trim_reads=fastqFiles
@jobs_limit(salmon_batch)
@follows(run_salmon_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/salmon')
@check_if_uptodate(check_file_exists)
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'salmon/{basename[0]}')
def run_salmon(input_files,output_folder):
    salmon(input_files,output_folder,salmonDb,thread,lib_type)
    if gff_type != '':
        print 'merge transcript expression to gene level'
        get_gene_expression(output_folder+'/quant.sf',gff_fn,gff_type)
#--------------------- 7. return finish message -----------------------------------------------------
@follows(run_salmon)
def last_function():
    Message('salmon finished',contact)

if __name__ == '__main__':
    try:
        pipeline_run([run_salmon],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
                 touch_files_only=False,verbose=5)
    except:
        Message('salmon failed',contact)
    