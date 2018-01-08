from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import hisat2_Db,hisat2
from Modules.Trimmomatic import conda_Trimmomatic
from Modules.Samtools import sortBam
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
rRNA_fa = p.rRNA_fa
ref_fa = p.ref_fa
annotation = p.gff
# trimmomatic parameter
trim = p.trim_reads
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter
# star parameter
hisat2_batch = p.hisat2_jobs_per_batch
hisat2_rRNA_db = p.hisat2_rrna_index
hisat2_target_db = p.hisat2_target_index

other = p.other
contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
Message('Riboseq start',contact)
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
    conda_Trimmomatic(input_file,output_file,n,adapter,20)
    remove(input_file)
    # run fastqc after trimming
    if QC:
        for fq in output_file:
            sarge.run('fastqc {input} -o fastqc'.format(input=fq))
#--------------------- 3. Align to rtRNA -----------------------------------------------
# build index
@active_if(not os.path.exists(hisat2_rRNA_db))
@follows(trim_reads,run_QC1)
def hisat2_rrna_index():
    if not os.path.exists(hisat2_rRNA_db): os.mkdir(hisat2_rRNA_db)
    hisat2_Db(rRNA_fa,hisat2_rRNA_db+'/rRNA',thread)
# align
if trim == False:
    trim_reads=fastqFiles
@jobs_limit(hisat2_batch)
@follows(hisat2_rrna_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/f01rRNA_bam')
@check_if_uptodate(check_file_exists)
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'f01rRNA_bam/{basename[0]}_norrna.fq.gz')
def run_hisat2rRNA(input_file,output_file):
    n = num_thread2use(hisat2_batch,len(fastqFiles),thread)
    rrna_fq = output_file[:-13]+'.bam'
    hisat2(input_file,rrna_fq,hisat2_rRNA_db+'/rRNA',n,['--un-gz',output_file])
#--------------------- 3. Align to target genome ----------------------------------------
@active_if(not os.path.exists(hisat2_target_db))
@follows(run_hisat2rRNA)
def hisat2_index():
    if not os.path.exists(hisat2_target_db): os.mkdir(hisat2_target_db)
    hisat2_Db(ref_fa,hisat2_target_db+'/target',thread)
# align
@jobs_limit(hisat2_batch)
@follows(hisat2_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/f02_bam')
@check_if_uptodate(check_file_exists)
@transform(run_hisat2rRNA,formatter('.*\.f.*?\.gz'),'f02_bam/{basename[0]}.bam')
def run_hisat2(input_file,output_file):
    n = num_thread2use(hisat2_batch,len(fastqFiles),thread)
    hisat2([input_file],output_file,hisat2_target_db+'/target',n,['--known-splicesite-infile',annotation])
#--------------------- 4. get primary mapping -------------------------------------------
@jobs_limit(trim_batch)
@follows(run_hisat2)
@mkdir(fastqFiles,formatter(),'{path[0]}/f03_primaryBam')
@check_if_uptodate(check_file_exists)
@transform(run_hisat2,formatter('.*\.bam'),'f03_primaryBam/{basename[0]}.bam')
def primary_bam(input_file,output_file):
    cmd = ('samtools view -h {fst_map} | grep -E {pattern} | '
               'samtools view -bh -F 256 - > {out}').format(fst_map=input_file,
                pattern='\'(NM:i:[012])|(^@)\'',out=output_file)
    print(cmd)
    sarge.run(cmd)
#--------------------- 5. samtools sort by position --------------------------------------
@jobs_limit(trim_batch)
@follows(primary_bam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f04_sortBam')
@check_if_uptodate(check_file_exists)
@transform(primary_bam,formatter('.*\.bam'),'f04_sortBam/{basename[0]}.sort.bam')
def sort_by_pos(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n)
    stat = sarge.get_stdout('samtools flagstat {bam}'.format(bam=output_file))
    with open(output_file[:-3]+'flagstat.txt','w') as f:
        f.write(stat)
# @follows(sort_by_pos)
# def remove_bam():
#     if os.path.exists('f01rRNA_bam'): shutil.rmtree('f01rRNA_bam')   # remove bam folder
 
@follows(sort_by_pos)
def last_function():
    Message('Riboseq finished',contact)

if __name__ == '__main__':
    try:
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
                 touch_files_only=False,verbose=15)
    except:
        Message('Riboseq failed',contact)
