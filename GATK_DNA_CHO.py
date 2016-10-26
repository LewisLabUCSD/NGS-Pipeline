from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import bwa_Db,bwa_mem
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import sortBam,build_fa_index,merge_bams
from Modules.Picard import build_fa_dict,mark_duplicates
from Modules.GATK import *
import yaml
import sys,shutil
import glob
from natsort import natsorted


#============ parameters ======================
parameter_file =  sys.argv[1]
#parameter_file = '/data/shangzhong/DE/winzeler/RNAseq_count.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)
#------------- get parameters -----------
file_path = p.RawDataPath
thread = p.thread
# all parameter
ref_fa = p.ref_fa
annotation = p.gff
# trimmomatic parameter
trim = p.trim_reads
trimmomatic = p.trimmomatic_path
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter

picard = p.picard
gatk = p.gatk

bwa_batch = p.bwa_jobs_per_batch
bwa_Db = p.bwa_Db

sp = p.sample_name
read_groups = p.read_groups

contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
#Message('GATK_DNA_CHO start',contact)
os.chdir(file_path)
#===============================================================================
#                     Part I. Preprocess
#===============================================================================
#--------------------- 1. build index for fa file using samtools and GATK ------------------
dict_file = '.'.join(ref_fa.split('.')[:-1]) + '.dict'
fai_file = ref_fa + '.fai'
if not os.path.exists(dict_file): build_fa_dict(ref_fa,picard)
if not os.path.exists(fai_file): build_fa_index(ref_fa)
#--------------------- 2. read all files ------------------------------------------------
fastqFiles = list_fq_files(file_path)
def trim_parameters():
    infiles,outfiles = replace_filename(fastqFiles,'^','trim_')
    for infile, output in zip(infiles,outfiles):
        yield infile,output
#--------------------- 3. trim reads-----------------------------------------------------
@active_if(trim)
@jobs_limit(trim_batch)
@files(trim_parameters)
def trim_reads(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
    remove(input_file)
#--------------------- 4. Map with bwa -----------------------------------------------------
def get_fq_and_readgroup():
    fqFiles = list_fq_files(file_path)
    for fq, rg in zip(fqFiles,read_groups):
        out = 'bam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
        yield fq,out,rg
# build index
@active_if(not os.path.exists(bwa_Db))
@follows(trim_reads)
def bwa_index():
    bwa_Db(bwa_Db,ref_fa)
# align
@jobs_limit(bwa_batch)
@follows(bwa_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/bam')
@files(get_fq_and_readgroup)
@check_if_uptodate(check_file_exists)
def run_bwa(input_file,output_file,rg):
    n = num_thread2use(bwa_batch,len(fastqFiles),thread)
    bwa_mem(input_file,output_file,bwa_Db+'/bwa',n,otherParameters=['-M','-R '+rg+'\\tPL:illumina\\tLB:lib20000\\tPU:unit1'])
#--------------------- 5. Sort bam file --------------------------------------------------
@jobs_limit(trim_batch)
@follows(run_bwa)
@mkdir(fastqFiles,formatter(),'{path[0]}/sortBam')
@transform(run_bwa,formatter('.*\.bam'),'sortBam/{basename[0]}.sort.bam')
@check_if_uptodate(check_file_exists)
def sort_by_pos(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='pos')
@follows(sort_by_pos)
def remove_bam():
    if os.path.exists('bam'): shutil.rmtree('bam')   # remove bam folder
#--------------------- 6. Markduplicates using picard -------------------------------------
@jobs_limit(trim_batch)
@follows(sort_by_pos)
@mkdir(fastqFiles,formatter(),'{path[0]}/dedupBam')
@transform(sort_by_pos,formatter('.*\.sort\.bam'),'dedupBam/{basename[0]}.dedup.bam')
@check_if_uptodate(check_file_exists)
def markduplicates(input_file,output_file):
    mark_duplicates(input_file,output_file,picard)
@follows(markduplicates)
def remove_sortBam():
    if os.path.exists('sortBam'): shutil.rmtree('sortBam')
#--------------------- 7. Indel realignment ---------------------
@jobs_limit(thread)
@follows(markduplicates)
@mkdir(fastqFiles,formatter(),'{path[0]}/indelReali')
@transform(markduplicates,formatter('.*\.dedup\.bam'),'indelReali/{basename[0]}.reali.bam')
@check_if_uptodate(check_file_exists)
def Realign(input_file,output_file):
    interval = re.sub('reali\.bam$','interval.list',output_file)
    RealignerTargetCreator(input_file,interval,gatk,ref_fa,1,gold_indels=[''])
    IndelRealigner(input_file,output_file,gatk,ref_fa,interval,gold_indels=[''])
@follows(Realign)
def remove_dedupBam():
    if os.path.exists('dedupBam'): shutil.rmtree('dedupBam')
#--------------------- 8. Round 1 call ----------------------------------------------
@jobs_limit(thread)
@follows(Realign)
@mkdir(fastqFiles,formatter(),'{path[0]}/Round1Call')
@transform(Realign,formatter('.*\.reali\.bam'),'Round1Call/{basename[0]}.raw.vcf')
@check_if_uptodate(check_file_exists)
def round1Vari_call(input_file,output_file):
    n = num_thread2use(thread,len(fastqFiles),thread)
    HaplotypeCaller_DNA_gVCF(input_file,output_file,gatk,ref_fa,n,otherParameters=[])
#--------------------- 9. Merge raw vcf ----------------------------------------------
@follows(round1Vari_call)
@merge(round1Vari_call,'Round1Call/round1.g.vcf')
@check_if_uptodate(check_file_exists)
def merge_vcf(input_files,output_file):
    JointGenotype(input_files,output_file,gatk,ref_fa,thread)
#--------------------- 10. filter gold snp and indel ---------------------------------
@follows(merge_vcf)
@transform(merge_vcf,suffix('.g.vcf'),['.gold_snp.vcf','.gold_indel.vcf'])
@check_if_uptodate(check_file_exists)
def hard_filter(input_file,output_pair):
    HardFilter(input_file,output_pair,gatk,ref_fa,thread)
#--------------------- 11. Base recalibration -----------------------------------------
# step 1
@follows(hard_filter)
@mkdir(fastqFiles,formatter(),'{path[0]}/BaseRecal')
@transform(Realign,formatter('.*\.reali\.bam'),add_inputs(hard_filter),'BaseRecal/{basename[0]}.table')
@check_if_uptodate(check_file_exists)
#@check_if_uptodate(check_file_exists)
def Baserecalibration_1(input_file,output_file):
    BaseRecalibrator_1(input_file[0],output_file,input_file[1],gatk,ref_fa,thread)
# step 2
@follows(Baserecalibration_1)
@transform(Realign,formatter('.*\.reali\.bam'),add_inputs(hard_filter),'BaseRecal/{basename[0]}.post_table')
@check_if_uptodate(check_file_exists)
def Baserecalibration_2(input_file,output_file):
    bam = input_file[0].split('/')[-1]
    table = 'BaseRecal/' + re.sub('\.reali\.bam$','table',bam)
    BaseRecalibrator_2(input_file[0],output_file,table,input_file[1],gatk,ref_fa,thread)
# step 3
@follows(Baserecalibration_2)
@transform(Baserecalibration_1,formatter('.+\.table'),'BaseRecal/{basename[0]}.plot')
@check_if_uptodate(check_file_exists)
def Baserecalibration_3(input_file,output_file):
    post_table = re.sub('\.table$','.post_table',input_file)
    BaseRecalibrator_3(input_file,output_file,post_table,gatk,ref_fa)
# step 4
@jobs_limit(trim_batch)
@follows(Baserecalibration_3)
@transform(Realign,formatter('.*\.reali\.bam'),add_inputs(hard_filter),'BaseRecal/{basename[0]}.recal.bam')
@check_if_uptodate(check_file_exists)
def Baserecalibration_4(input_file,output_file):
    table = re.sub('\.bam','.table',input_file[0])
    BaseRecalibrator_4(input_file[0],output_file,gatk,ref_fa,input_file[1],table,thread)
@follows(Baserecalibration_4)
def remove_realiBam():
    if os.path.exists('indelReali'): shutil.rmtree('indelReali')
#--------------------- 12. merge lanes for the same sample -----------------------------------------
def get_group_bam():
    readic = {}
    bamfiles = natsorted(glob.glob('BaseRecal/*.recal.bam'))
    for rg,bam in zip(read_groups,bamfiles):
        start = rg.index('SM:')
        sample = rg[start+3:]
        if sample in readic:
            readic[sample].append(bam)
        else:
            readic[sample] = [bam]
    for sp in readic:
        output_file = 'mergeBam/' + sp + '.merge.bam'
        input_file = readic[sp]
        yield input_file,output_file
@follows(Baserecalibration_4)
@mkdir(fastqFiles,formatter(),'{path[0]}/mergeBam')
@files(get_group_bam)
@check_if_uptodate(check_file_exists)
def mergeBam(input_files,output_file):
    merge_bams(input_files,output_file)
#--------------------- 13. Mark duplicates for merged file ---------------------------------------
@jobs_limit(trim_batch)
@follows(mergeBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/dedupBam2')
@transform(mergeBam,formatter('.*\.merge\.bam'),'dedupBam2/{basename[0]}.dedup.bam')
@check_if_uptodate(check_file_exists)
def markduplicates2(input_file,output_file):
    mark_duplicates(input_file,output_file,picard)
@follows(markduplicates2)
def remove_mergeBam():
    if os.path.exists('mergeBam'): shutil.rmtree('mergeBam')
#--------------------- 14. Indel realignment -----------------------------------------------------
@jobs_limit(thread)
@follows(markduplicates2)
@mkdir(fastqFiles,formatter(),'{path[0]}/indelReali2')
@transform(markduplicates2,formatter('.*\.dedup\.bam'),'indelReali2/{basename[0]}.reali.bam')
@check_if_uptodate(check_file_exists)
def Realign2(input_file,output_file):
    interval = re.sub('reali\.bam$','interval.list',output_file)
    RealignerTargetCreator(input_file,interval,gatk,ref_fa,1,gold_indels=[''])
    IndelRealigner(input_file,output_file,gatk,ref_fa,interval,gold_indels=[''])
@follows(Realign2)
def remove_dedupBam2():
    if os.path.exists('dedupBam2'): shutil.rmtree('dedupBam2')
#--------------------- 15. Round 2 call ----------------------------------------------
@follows(Realign2)
@mkdir(fastqFiles,formatter(),'{path[0]}/Round2Call')
@transform(Realign2,formatter('.*\.reali\.bam'),'Round2Call/{basename[0]}.raw.vcf')
@check_if_uptodate(check_file_exists)
def round2Vari_call(input_file,output_file):
    n = num_thread2use(thread,len(fastqFiles),thread)
    HaplotypeCaller_DNA_gVCF(input_file,output_file,gatk,ref_fa,n,otherParameters=[])
#--------------------- 16. Merge raw2 vcf ---------------------------------------------
@follows(round2Vari_call)
@merge(round2Vari_call,'Round2Call/round2.g.vcf')
@check_if_uptodate(check_file_exists)
def merge_vcf2(input_files,output_file):
    JointGenotype(input_files,output_file,gatk,ref_fa,thread)
#--------------------- 17. filter gold snp and indel ---------------------------------
@follows(merge_vcf2)
@transform(merge_vcf2,suffix('.g.vcf'),['.gold_snp.vcf','.gold_indel.vcf'])
@check_if_uptodate(check_file_exists)
def hard_filter2(input_file,output_pair):
    HardFilter(input_file,output_pair,gatk,ref_fa,thread)
#--------------------- 18. combine vcf files ---------------------------------
@follows(hard_filter2)
@merge(hard_filter2,sp+'.vcf')
def combine_vcf(input_files,output_file):
    CombineSNPandINDEL(input_files,output_file,gatk,ref_fa,otherParams=['--assumeIdenticalSamples','--genotypemergeoption UNSORTED'])

@follows(combine_vcf)
def last_function():
    Message('test succeed',contact)
    

if __name__ == '__main__':
    try:
#         pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
        touch_files_only=False,verbose=5)
    except:
        Message('test failed',contact)

    
    
    
    
    
    
    
    