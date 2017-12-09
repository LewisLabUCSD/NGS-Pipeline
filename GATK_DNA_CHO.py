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
#parameter_file = '/data/shangzhong/DNArepair/fq/GATK_DNA_CHO.yaml'
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

QC = p.QC
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
Message('GATK_DNA_start',contact)
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
if fastqFiles[0][0].startswith('trim_'):
    trim = False
def trim_parameters():
    infiles,outfiles = replace_filename(fastqFiles,'^','trim_')
    for infile, output in zip(infiles,outfiles):
        yield infile,output
#--------------------- run fastqc before trimming -----------
@active_if(QC)
@jobs_limit(thread)
@mkdir(fastqFiles,formatter(),'{path[0]}/fastqc')
@files(trim_parameters)
def run_QC1(input_file,output_file):
    for fq in input_file:
        sarge.run('fastqc {input} -o fastqc'.format(input=fq))
#---------------------3. trim file ------------------
@active_if(trim)
@follows(run_QC1)
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
        out = 'f01_bam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
        yield fq,out,rg
# build index
@active_if(not os.path.exists(bwa_Db))
@follows(trim_reads)
def bwa_index():
    bwa_Db(bwa_Db,ref_fa)
# align
@jobs_limit(bwa_batch)
@follows(bwa_index,trim_reads,run_QC1)
@mkdir(fastqFiles,formatter(),'{path[0]}/f01_bam')
@files(get_fq_and_readgroup)
def run_bwa(input_file,output_file,rg):
    n = num_thread2use(bwa_batch,len(fastqFiles),thread)
    bwa_mem(input_file,output_file,bwa_Db+'/bwa',n,otherParameters=['-R '+rg+'\\\\tPL:illumina\\\\tLB:lib20000\\\\tPU:unit1'])
#--------------------- 5. Sort bam file --------------------------------------------------
@jobs_limit(trim_batch)
@follows(run_bwa)
@mkdir(fastqFiles,formatter(),'{path[0]}/f02_sortBam')
@transform(run_bwa,formatter('.*\.bam'),'f02_sortBam/{basename[0]}.sort.bam')
@check_if_uptodate(check_file_exists)
def sort_by_pos(input_file,output_file):
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    sortBam(input_file,output_file,n,sortType='pos')
@follows(sort_by_pos)
def remove_bam():
    if os.path.exists('f01_bam'): shutil.rmtree('f01_bam')   # remove bam folder
#--------------------- 6. Markduplicates using picard -------------------------------------
@jobs_limit(trim_batch)
@follows(remove_bam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f03_dedupBam')
@transform(sort_by_pos,formatter('.*\.sort\.bam'),'f03_dedupBam/{basename[0]}.dedup.bam')
@check_if_uptodate(check_file_exists)
def markduplicates(input_file,output_file):
    mark_duplicates(input_file,output_file,picard)
@follows(markduplicates)
def remove_sortBam():
    if os.path.exists('f02_sortBam'): shutil.rmtree('f02_sortBam')
#--------------------- 7. Indel realignment ---------------------
@jobs_limit(thread)
@follows(remove_sortBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f04_indelReali')
@transform(markduplicates,formatter('.*\.dedup\.bam'),'f04_indelReali/{basename[0]}.reali.bam')
@check_if_uptodate(check_file_exists)
def Realign(input_file,output_file):
    interval = re.sub('reali\.bam$','interval.list',output_file)
    RealignerTargetCreator(input_file,interval,gatk,ref_fa,1,gold_indels=[''])
    IndelRealigner(input_file,output_file,gatk,ref_fa,interval,gold_indels=[''])
@follows(Realign)
def remove_dedupBam():
    if os.path.exists('f03_dedupBam'): shutil.rmtree('f03_dedupBam')
#--------------------- 8. Round 1 call ----------------------------------------------
@jobs_limit(thread)
@follows(remove_dedupBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f05_Round1Call')
@transform(Realign,formatter('.*\.reali\.bam'),'f05_Round1Call/{basename[0]}.raw.g.vcf')
@check_if_uptodate(check_file_exists)
def round1Vari_call(input_file,output_file):
    n = num_thread2use(thread,len(fastqFiles),thread)
    HaplotypeCaller_DNA_gVCF(input_file,output_file,gatk,ref_fa,n,otherParameters=[])
#--------------------- 9. Merge raw vcf ----------------------------------------------
@follows(round1Vari_call)
@merge(round1Vari_call,'f05_Round1Call/round1.g.vcf')
@check_if_uptodate(check_file_exists)
def merge_vcf(input_files,output_file):
    JointGenotype(input_files,output_file,gatk,ref_fa,thread)
#--------------------- 10. filter gold snp and indel ---------------------------------
@follows(merge_vcf)
@transform(merge_vcf,suffix('.g.vcf'),['.gold_snp.vcf','.gold_indel.vcf'])
# @check_if_uptodate(check_file_exists)
def hard_filter(input_file,output_pair):
    HardFilter(input_file,output_pair,gatk,ref_fa,thread)
@follows(hard_filter)
def remove_vcf():
    if os.path.exists('f05_Round1Call'):
        for f in glob.glob('f05_Round1Call/*'):
            if 'gold' not in f:
                os.remove(f)
#--------------------- 11. Base recalibration -----------------------------------------
# step 1
@follows(remove_vcf)
@mkdir(fastqFiles,formatter(),'{path[0]}/f06_BaseRecal')
@transform(Realign,formatter('.*\.reali\.bam'),add_inputs(hard_filter),'f06_BaseRecal/{basename[0]}.table')
@check_if_uptodate(check_file_exists)
def Baserecalibration_1(input_file,output_file):
    n = num_thread2use(thread,len(fastqFiles),thread)
    BaseRecalibrator_1(input_file[0],output_file,input_file[1],gatk,ref_fa,thread=str(n))
# step 2
@follows(Baserecalibration_1)
@transform(Realign,formatter('.*\.reali\.bam'),add_inputs(hard_filter),'f06_BaseRecal/{basename[0]}.post_table')
@check_if_uptodate(check_file_exists)
def Baserecalibration_2(input_file,output_file):
    bam = input_file[0].split('/')[-1]
    table = 'f06_BaseRecal/' + re.sub('\.reali\.bam$','.reali.table',bam)
    n = num_thread2use(thread,len(fastqFiles),thread)
    BaseRecalibrator_2(input_file[0],output_file,table,input_file[1],gatk,ref_fa,thread=str(n))
# step 3
@follows(Baserecalibration_2)
@transform(Baserecalibration_1,formatter('.+\.table'),'f06_BaseRecal/{basename[0]}.plot')
@check_if_uptodate(check_file_exists)
def Baserecalibration_3(input_file,output_file):
    post_table = re.sub('\.table$','.post_table',input_file)
    BaseRecalibrator_3(input_file,output_file,post_table,gatk,ref_fa)
# step 4
@jobs_limit(bwa_batch)
@follows(Baserecalibration_3)
@transform(Realign,formatter('.*\.reali\.bam'),add_inputs(hard_filter),'f06_BaseRecal/{basename[0]}.recal.bam')
@check_if_uptodate(check_file_exists)
def Baserecalibration_4(input_file,output_file):
    table = 'f06_BaseRecal/'+re.sub('\.bam','.table',input_file[0].split('/')[1])
    n = num_thread2use(bwa_batch,len(fastqFiles),thread)
    BaseRecalibrator_4(input_file[0],output_file,gatk,ref_fa,input_file[1],table,n)
@follows(Baserecalibration_4)
def remove_realiBam():
    if os.path.exists('f04_indelReali'): shutil.rmtree('f04_indelReali')
#--------------------- 12. merge lanes for the same sample -----------------------------------------
def get_rg_dic():
    readic = {}
    bamfiles = natsorted(glob.glob('f06_BaseRecal/*.recal.bam'))
    for rg,bam in zip(read_groups,bamfiles):
        start = rg.index('SM:')
        sample = rg[start+3:]
        if sample in readic:
            readic[sample].append(bam)
        else:
            readic[sample] = [bam]
    return readic
def get_group_bam():
    readic = get_rg_dic()
    for sp in readic:
        output_file = 'f07_mergeBam/' + sp + '.merge.bam'
        input_file = readic[sp]
        yield input_file,output_file
@follows(remove_realiBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f07_mergeBam')
@files(get_group_bam)
def mergeBam(input_files,output_file):
    merge_bams(input_files,output_file)
@follows(mergeBam)
def remove_recalBam():
    if os.path.exists('f06_BaseRecal'):
        for f in glob.glob('f06_BaseRecal/*'):
            os.remove(f)
            if f.endswith('recal.bam'):
                handle = open(f,'w')
                handle.close()
#--------------------- 13. Mark duplicates for merged file ---------------------------------------
@jobs_limit(trim_batch)
@follows(remove_recalBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f08_dedupBam2')
@transform(mergeBam,formatter('.*\.merge\.bam'),'f08_dedupBam2/{basename[0]}.dedup.bam')
@check_if_uptodate(check_file_exists)
def markduplicates2(input_file,output_file):
    mark_duplicates(input_file,output_file,picard)
@follows(markduplicates2)
def remove_mergeBam():
    if os.path.exists('f07_mergeBam'): shutil.rmtree('f07_mergeBam')
#--------------------- 14. Indel realignment -----------------------------------------------------
@follows(remove_mergeBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f09_indelReali2')
@transform(markduplicates2,formatter('.*\.dedup\.bam'),'f09_indelReali2/{basename[0]}.reali.bam')
@check_if_uptodate(check_file_exists)
def Realign2(input_file,output_file):
    interval = re.sub('reali\.bam$','interval.list',output_file)
    n = num_thread2use(len(get_rg_dic().keys()),len(fastqFiles),thread)
    RealignerTargetCreator(input_file,interval,gatk,ref_fa,n,gold_indels=[''])
    IndelRealigner(input_file,output_file,gatk,ref_fa,interval,gold_indels=[''])
@follows(Realign2)
def remove_dedupBam2():
    if os.path.exists('f08_dedupBam2'): shutil.rmtree('f08_dedupBam2')
#--------------------- 15. Round 2 call ----------------------------------------------
@follows(remove_dedupBam2)
@mkdir(fastqFiles,formatter(),'{path[0]}/f10_Round2Call')
@transform(Realign2,formatter('.*\.reali\.bam'),'f10_Round2Call/{basename[0]}.raw.g.vcf')
@check_if_uptodate(check_file_exists)
def round2Vari_call(input_file,output_file):
    n = num_thread2use(len(get_rg_dic().keys()),len(fastqFiles),thread)
    HaplotypeCaller_DNA_gVCF(input_file,output_file,gatk,ref_fa,n,otherParameters=[])
#--------------------- 16. Merge raw2 vcf ---------------------------------------------
@follows(round2Vari_call)
@merge(round2Vari_call,'f10_Round2Call/round2.g.vcf')
def merge_vcf2(input_files,output_file):
    JointGenotype(input_files,output_file,gatk,ref_fa,thread)
#--------------------- 17. filter gold snp and indel ---------------------------------
@follows(merge_vcf2)
@transform(merge_vcf2,suffix('.g.vcf'),['.gold_snp.vcf','.gold_indel.vcf'])
# @check_if_uptodate(check_file_exists)
def hard_filter2(input_file,output_pair):
    HardFilter(input_file,output_pair,gatk,ref_fa,thread)
#--------------------- 18. combine vcf files ---------------------------------
@follows(hard_filter2)
@mkdir(fastqFiles,formatter(),'{path[0]}/f11_FinalVcf')
@merge(hard_filter2,'f11_FinalVcf/'+sp+'.merged.filter.vcf')
@check_if_uptodate(check_file_exists)
def combine_vcf(input_files,output_file):
    CombineSNPandINDEL(input_files,output_file,gatk,ref_fa,otherParams=['--assumeIdenticalSamples','--genotypemergeoption UNSORTED'])

@follows(combine_vcf)
def last_function():
    Message('GATK_DNA_succeed',contact)
    

if __name__ == '__main__':
    try:
#         pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
        touch_files_only=False,verbose=5)
    except:
        pass
        Message('test failed',contact)

    
    


    
    
    
    
    