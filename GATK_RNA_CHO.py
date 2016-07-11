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
#parameter_file = '/data/shangzhong/Proteogenomics/test/GATK_RNA_CHO.yaml'
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

star_batch = p.star_jobs_per_batch
star_db = p.star_index

sp = p.sample_name
read_groups = p.read_groups

contact = p.contact
#===============================================================================
#                    Pipeline part
#===============================================================================
Message('GATK_DNA_CHO start',contact)
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
#--------------------- 4. Map with STAR -----------------------------------------------------
# build index
@active_if(not os.path.exists(star_db))
@follows(trim_reads)
def star_index():
    STAR_Db(star_db,ref_fa,thread)
# align
if trim == True:
    @jobs_limit(star_batch)
    @follows(trim_reads,star_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/f01_bam')
    @check_if_uptodate(check_file_exists)
    @transform(trim_reads,formatter('.*\.f.*?\.gz'),'f01_bam/{basename[0]}.sort.bam')
    def run_star(input_file,output_file):
        n = num_thread2use(star_batch,len(fastqFiles),thread)
        STAR(input_file,output_file,star_db,n,'',['--outSAMtype BAM','SortedByCoordinate','--twopassMode Basic'])
else:
    @jobs_limit(star_batch)
    @follows(star_index)
    @mkdir(fastqFiles,formatter(),'{path[0]}/f01_bam')
    @check_if_uptodate(check_file_exists)
    @transform(fastqFiles,formatter('.*\.f.*?\.gz'),'f01_bam/{basename[0]}.sort.bam')
    def run_star(input_file,output_file):
        n = num_thread2use(star_batch,len(fastqFiles),thread)
        STAR(input_file,output_file,star_db,n,'',['--outSAMtype BAM','SortedByCoordinate','--twopassMode Basic'])
#--------------------- 5. add read group --------------------------------------------------
def get_bam_and_rg():
    bams = [f for f in os.listdir('f01_bam') if f.endswith('.sort.bam')]
    bams = natsorted(bams)
    for bam, rg in zip(bams,read_groups):
        output = re.sub('\.sort\.bam','.adrg.bam',bam)
        yield ['f01_bam/'+bam,rg],'f02_addGroup/' + output
@jobs_limit(trim_batch*2)
@follows(run_star,mkdir('tmp'),mkdir('f02_addGroup'))
@files(get_bam_and_rg)
@check_if_uptodate(check_file_exists)
def run_add_group(input_file,output_file):
    add_readgroup(input_file[0],output_file,input_file[1],picard)
@follows(run_add_group)
def remove_bam():
    if os.path.exists('f01_bam'):shutil.rmtree('f01_bam')
#--------------------- 6. Markduplicates using picard -------------------------------------
@jobs_limit(trim_batch)
@follows(run_add_group,remove_bam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f03_dedupBam')
@transform(run_add_group,formatter('.*\.adrg\.bam'),'f03_dedupBam/{basename[0]}.dedup.bam')
@check_if_uptodate(check_file_exists)
def run_mark_duplicates(input_file,output_file):
    mark_duplicates(input_file,output_file,picard)
@follows(run_mark_duplicates)
def remove_groupBam():
    if os.path.exists('f02_addGroup'): shutil.rmtree('f02_addGroup')
#--------------------- 7. Split N ---------------------------------------------------------
@follows(run_mark_duplicates,remove_groupBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f04_splitBam')
@check_if_uptodate(check_file_exists)
@transform(run_mark_duplicates,formatter('.*\.dedup\.bam'),'f04_splitBam/{basename[0]}.split.bam')
def run_splitN(input_file,output_file):
    splitN(input_file,output_file,gatk,ref_fa)
@follows(run_splitN)
def remove_dedupBam():
    if os.path.exists('f03_dedupBam'): shutil.rmtree('f03_dedupBam')
#--------------------- 8. Indel realignment ---------------------
@jobs_limit(thread)
@follows(run_splitN,remove_dedupBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f05_indelReali')
@transform(run_splitN,formatter('.*\.split\.bam'),'f05_indelReali/{basename[0]}.reali.bam')
@check_if_uptodate(check_file_exists)
def run_realign(input_file,output_file):
    interval = re.sub('reali\.bam$','interval.list',output_file)
    RealignerTargetCreator(input_file,interval,gatk,ref_fa,1,gold_indels=[''])
    IndelRealigner(input_file,output_file,gatk,ref_fa,interval,gold_indels=[''])
@follows(run_realign)
def remove_splitBam():
    if os.path.exists('f04_splitBam'): shutil.rmtree('f04_splitBam')
#--------------------- 9. Round 1 call ----------------------------------------------
@jobs_limit(thread)
@follows(run_realign,remove_splitBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f06_Round1Call')
@transform(run_realign,formatter('.*\.reali\.bam'),'f06_Round1Call/{basename[0]}.vcf')
@check_if_uptodate(check_file_exists)
def round1Vari_call(input_file,output_file):
    n = num_thread2use(thread,len(fastqFiles),thread)
    HaplotypeCaller_RNA_VCF(input_file,output_file,gatk,ref_fa,n)
#--------------------- 10. filter gold snp and indel ---------------------------------
@follows(round1Vari_call)
@transform(round1Vari_call,suffix('.vcf'),'.gold.vcf')
@check_if_uptodate(check_file_exists)
def run_RNA_Vari_Filter(input_file,output_file):
    RNA_Vari_Filter(input_file,output_file,gatk,ref_fa)
#--------------------- 11. Base recalibration -----------------------------------------
# step 1
@follows(run_RNA_Vari_Filter)
@mkdir(fastqFiles,formatter(),'{path[0]}/f07_BaseRecal')
@transform(run_realign,formatter('.*\.reali\.bam'),'f07_BaseRecal/{basename[0]}.table')
@check_if_uptodate(check_file_exists)
def run_RNA_Baserecalibration_1(input_file,output_file):
    gold_vcf = 'f06_Round1Call/'+re.sub('\.bam$','.gold.vcf',input_file).split('/')[-1]
    RNA_BaseRecalibrator_1(input_file,output_file,gatk,ref_fa,gold_vcf,thread)
# step 2
@follows(run_RNA_Baserecalibration_1)
@transform(run_realign,formatter('.*\.reali\.bam'),'f07_BaseRecal/{basename[0]}.post_table')
@check_if_uptodate(check_file_exists)
def run_RNA_Baserecalibration_2(input_file,output_file):
    table = re.sub('post_table$','table',output_file)
    gold_vcf = 'f06_Round1Call/'+re.sub('\.bam$','.gold.vcf',input_file).split('/')[-1]
    RNA_BaseRecalibrator_2(input_file,output_file,table,gatk,ref_fa,gold_vcf,thread='1')
# step 3
@follows(run_RNA_Baserecalibration_2)
@transform(run_RNA_Baserecalibration_1,formatter('.+\.table'),'f07_BaseRecal/{basename[0]}.plot.pdf')
@check_if_uptodate(check_file_exists)
def run_RNA_Baserecalibration_3(input_file,output_file):
    post_table = re.sub('\.table$','.post_table',input_file)
    RNA_BaseRecalibrator3(input_file,output_file,post_table,gatk,ref_fa)
# step 4
@jobs_limit(trim_batch)
@follows(run_RNA_Baserecalibration_3)
@transform(run_realign,formatter('.*\.reali\.bam'),'f07_BaseRecal/{basename[0]}.recal.bam')
@check_if_uptodate(check_file_exists)
def run_RNA_Baserecalibration_4(input_file,output_file):
    table = re.sub('\.recal\.bam','.table',output_file)
    gold_vcf = 'f06_Round1Call/'+re.sub('\.bam$','.gold.vcf',input_file).split('/')[-1]
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    RNA_BaseRecalibrator4(input_file,output_file,gatk,table,ref_fa,gold_vcf,n)
@follows(run_RNA_Baserecalibration_4)
def remove_realiBam():
    if os.path.exists('f05_indelReali'): shutil.rmtree('f05_indelReali')
#--------------------- 12. merge lanes for the same sample -----------------------------------------
def get_group_bam():
    readic = {}
    bamfiles = natsorted(glob.glob('f07_BaseRecal/*.recal.bam'))
    for rg,bam in zip(read_groups,bamfiles):
        start = rg.index('SM:')
        sample = rg[start+3:]
        if sample in readic:
            readic[sample].append(bam)
        else:
            readic[sample] = [bam]
    for sp in readic:
        output_file = 'f08_mergeBam/' + sp + '.merge.bam'
        input_file = readic[sp]
        yield input_file,output_file
@follows(run_RNA_Baserecalibration_4,remove_realiBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f08_mergeBam')
@files(get_group_bam)
@check_if_uptodate(check_file_exists)
def mergeBam(input_files,output_file):
    merge_bams(input_files,output_file)
#--------------------- 13. Mark duplicates for merged file ---------------------------------------
@jobs_limit(trim_batch)
@follows(mergeBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f09_dedupBam2')
@transform(mergeBam,formatter('.*\.merge\.bam'),'f09_dedupBam2/{basename[0]}.dedup.bam')
@check_if_uptodate(check_file_exists)
def markduplicates2(input_file,output_file):
    mark_duplicates(input_file,output_file,picard)
@follows(markduplicates2)
def remove_mergeBam():
    if os.path.exists('f08_mergeBam'): shutil.rmtree('f08_mergeBam')
#--------------------- 14. Indel realignment -----------------------------------------------------
@jobs_limit(thread)
@follows(markduplicates2,remove_mergeBam)
@mkdir(fastqFiles,formatter(),'{path[0]}/f10_indelReali2')
@transform(markduplicates2,formatter('.*\.dedup\.bam'),'f10_indelReali2/{basename[0]}.reali.bam')
@check_if_uptodate(check_file_exists)
def run_realign2(input_file,output_file):
    interval = re.sub('reali\.bam$','interval.list',output_file)
    RealignerTargetCreator(input_file,interval,gatk,ref_fa,1,gold_indels=[''])
    IndelRealigner(input_file,output_file,gatk,ref_fa,interval,gold_indels=[''])
@follows(run_realign2)
def remove_dedupBam2():
    if os.path.exists('f09_dedupBam2'): shutil.rmtree('f09_dedupBam2')
#--------------------- 15. Round 2 call ----------------------------------------------
@follows(run_realign2,remove_dedupBam2)
@mkdir(fastqFiles,formatter(),'{path[0]}/f11_Round2Call')
@transform(run_realign2,formatter('.*\.reali\.bam'),'f11_Round2Call/{basename[0]}.vcf')
@check_if_uptodate(check_file_exists)
def run_round2Vari_call(input_file,output_file):
    n = num_thread2use(thread,len(fastqFiles),thread)
    HaplotypeCaller_RNA_VCF(input_file,output_file,gatk,ref_fa,n)
#--------------------- 16. filter final vcf ---------------------------------
@jobs_limit(thread)
@follows(run_round2Vari_call)
@mkdir(fastqFiles,formatter(),'{path[0]}/f12_FinalVcf')
@transform(run_round2Vari_call,formatter('.*\.vcf'), 'f12_FinalVcf/{basename[0]}.merged.filter.vcf')
@check_if_uptodate(check_file_exists)
def run_filter_2(input_file,output_file):
    output_file = output_file.split('.')[0] + '.merged.filter.vcf'
    RNA_Vari_Filter(input_file,output_file,gatk,ref_fa)
    
@follows(run_filter_2)
def last_function():
    Message('test succeed',contact)
    
if __name__ == '__main__':
    try:
#         pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True, 
        touch_files_only=False)
    except:
        Message('test failed',contact)

    
    
    
    
    
    
    
    