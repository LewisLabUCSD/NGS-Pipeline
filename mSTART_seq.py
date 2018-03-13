from ruffus import *
from Modules.f01_file_process import *
from Modules.Trimmomatic import Trimmomatic
from Modules.Aligner import STAR,STAR_Db
from Modules.Samtools import sortBam
from Modules.Homer import *
import yaml
import shutil
import itertools

#============ parameters ======================
parameter_file =  sys.argv[1]
# parameter_file = '/data/isshamie/TSS_CHO/mSTART'
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
@mkdir(fastqFiles,formatter(),'{path[0]}/f01_fastqc')
@files(trim_parameters)
def run_QC1(input_file,output_file):
    for fq in input_file:
        sarge.run('fastqc {input} -o f01_fastqc'.format(input=fq))
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
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'f01_fastqc/{basename[0]}')
@check_if_uptodate(check_file_exists)
# @files(trim_parameters)
def run_QC2(input_file,output_file):
    for fq_in,fq_out in zip(input_file,output_file):
        if fq_in.startswith('trim_'):
            sarge.run('fastqc {input} -o f01_fastqc'.format(input=fq_in))
        else:
            sarge.run('fastqc {input} -o f01_fastqc'.format(input=fq_out))
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
@mkdir(fastqFiles,formatter(),'{path[0]}/f02_bam')
@mkdir(fastqFiles,formatter(),'{path[0]}/f02_flagstat')
@check_if_uptodate(check_file_exists)
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'f02_bam/{basename[0]}.bam')
def run_star(input_file,output_file):
    n = num_thread2use(star_batch,len(fastqFiles),thread)
    STAR(input_file,output_file,db_path,n,annotation,['--outSAMtype BAM','SortedByCoordinate','--outSAMunmapped Within'])
    stat = sarge.get_stdout('samtools flagstat {bam}'.format(bam=output_file))
    flag_fn = output_file[:-3]+'flagstat.txt'
    with open(flag_fn,'w') as f:
        f.write(stat)
    shutil.move(flag_fn,'f02_flagstat')
#--------------------- 4. make tag_directory ------------------------------------------------------
@follows(run_star)
@mkdir(fastqFiles,formatter(),'{path[0]}/f03_tags')
@check_if_uptodate(check_file_exists)
@transform(run_star,formatter('\.bam'),'f03_tags/{basename[0]}')
def make_tag(input_bam,out_dir):
    make_tag_directory(input_bam,out_dir,ref_fa)
    hist_out = out_dir+'/hist.txt'
    hist(out_dir,hist_out,ref_fa,annotation,mode='tss',peak='',region=4000,res=10,pc=3)
    hist_plot(hist_out)
#--------------------- 5. find peaks ------------------------------------------------------
def get_input_for_peak_call():
    gro_cap = [f for f in os.listdir('f03_tags') if 'mSTART' in f and 'input' not in f]
#     gro_cap_ctr = [[f for f in os.listdir('f03_tags') if 'contr' in f]]
    gro_seq = [f for f in os.listdir('f03_tags') if 'input' in f]
    comb = list(itertools.product(gro_cap,gro_seq))
    for com in comb:
        out = com[0] + '_and_' + com[1]
        yield ['f03_tags/'+f for f in com],'f04_peaks/' + out + '.peak'
@jobs_limit(thread)
@follows(make_tag)
@mkdir(fastqFiles,formatter(),'{path[0]}/f04_peaks')
@files(get_input_for_peak_call)
def find_peak(input_files,output_file):
    find_peaks(input_files[0],output_file,'tss',input_files[1],['-F 2'])
# #--------------------- 6. merge peaks ------------------------------------------------------
# @follows(find_peak)
# @merge(find_peak,'f04_peaks/merge.peak')
# def merge_peak(input_files,output_file):
#     merge_peaks(input_files,output_file,150)
# #--------------------- 6. annotate peaks ------------------------------------------------------
# @jobs_limit(thread)
# @follows(merge_peak)
# @mkdir(fastqFiles,formatter(),'{path[0]}/f05_annoPeaks')
# @transform(merge_peak,formatter('\.peak'),'f05_annoPeaks/{basename[0]}.anno')
# @check_if_uptodate(check_file_exists)
# def anno_peak(input_file,output_file):
#     annotate_peaks(input_file,output_file,ref_fa,annotation)
# #--------------------- 7. hist peaks ------------------------------------------------------
# @jobs_limit(thread)
# @follows(anno_peak)
# @mkdir(fastqFiles,formatter(),'{path[0]}/f06_histPeaks')
# @transform(find_peak,formatter('\.peak'),'f06_histPeaks/{basename[0]}.hist')
# @check_if_uptodate(check_file_exists)
# def peak_cov_hist(input_file,output_file): # input is peak file
#     gro_cap = [f for f in os.listdir('f03_tags') if '5GRO' in f]
#     tag = ['f03_tags/' + t for t in gro_cap if t in input_file]
#     hist(tag[0],output_file,ref_fa,annotation,mode='peak',peak=input_file,region=4000,res=25,pc=1)
#     hist_plot(output_file)

#--------------------- 8. merge peaks ------------------------------------------------------

# @follows(peak_cov_hist)
@follows(find_peak)
def last_function():
    Message('mSTART finished',contact)    


if __name__ == '__main__':
    try:
#         pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = False, 
                    touch_files_only=False,verbose=20)
    except:
        Message('mSTART failed',contact)
    

