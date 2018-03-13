from ruffus import *
from Modules.f01_file_process import *
from Modules.Aligner import STAR_Db,STAR
from Modules.Trimmomatic import Trimmomatic
from Modules.Samtools import *
from Modules.Homer import *
from Modules.peak_descriptives import *
import itertools
import yaml,sys
import shutil
import glob
import pickle
import inspect

#============ Load parameters ======================
parameter_file =  sys.argv[1]

#parameter_file = 'START_seq.yaml'
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = dic2obj(**doc)

#------------- get parameters -----------
file_path = p.RawDataPath #Folder contains the raw data and will create all the folders in there
thread = p.thread #Number of parallel threads to use 
contact = p.contact
QC = p.QC
# annotation parameters
ref_fa = p.ref_fa
annotation = p.annotation
# trimmomatic parameters
trim = p.trim_reads
trimmomatic = p.trimmomatic_path
is_homer_trim = p.is_homer_trim
trim_batch = p.trim_jobs_per_batch
adapter = p.adapter
# STAR parameters
star_batch = p.star_jobs_per_batch
star_db = p.star_index_path
run_pass = p.star_pass
other_params = p.star_params

# Homer parameters
peak_bg = p.homer_peak_background_fold #Need a background dataset to prevent False Positive peak aclls
type_merge = p.type_merge ##'all','one' how Homer merges replicates. All means they all have to have peaks

# Downstream parameters
multiple_peaks_width = p.multiple_peaks_width
divergent_width = p.divergent_width
# Downstream annotation parameters
mRNA_peak_file = p.mRNA_peak_file
marks = ['promoter','Intergenic']


#===============================================================================
#                    Pipeline part
#===============================================================================
Message('START-seq analysis starting',contact)
os.chdir(file_path)

##Save the parameters
with open('./params.txt','w') as f:
    f.write('\n'.join([i[0] + ':' + str(i[1]) for i in inspect.getmembers(p,lambda a:not(inspect.isroutine(a)))]))
 


#===============================================================================
#                     Part I. Preprocess
#===============================================================================
#--------------------- 1. read all files ------------------------------------------------
fastqFiles = list_fq_files(file_path)
print "Files we're working with:", fastqFiles

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
    print input_file,output_file
    n = num_thread2use(trim_batch,len(fastqFiles),thread)
    if is_homer_trim:
        homer_trim_cmd(input_file[0])
    else:
       Trimmomatic(input_file,output_file,trimmomatic,n,adapter)
    #remove(input_file)



#------------ run fastqc after trimming ------------
@active_if(QC and trim)
@follows(trim_reads)
@jobs_limit(thread)
@mkdir(fastqFiles,formatter(),'{path[0]}/f01_fastqc')
@transform(trim_reads,formatter('.*\.f.*?\.gz'),'f01_fastqc/{basename[0]}')
#@check_if_uptodate(check_file_exists)
#@files(trim_parameters)
def run_QC2(input_file,output_file):
    for fq_in,fq_out in zip(input_file,output_file):
        if fq_in.startswith('trim_'):
            sarge.run('fastqc {input} -o f01_fastqc'.format(input=fq_in))
        else:
            sarge.run('fastqc {input} -o f01_fastqc'.format(input=fq_out))


#--------------------- 3. Map with STAR -----------------------------------------------------
def get_fq():
    fqFiles = list_fq_files(file_path)
    for fq in fqFiles:
        out = 'f02_bam/' + re.sub('\.f.*q\.gz','.bam',fq[0])
        yield fq,out

# # build index
@active_if(not os.path.exists(star_db))
@follows(run_QC2)
def star_index():
    STAR_Db(star_db,ref_fa,thread)

# align
other_params.extend(['--outSAMtype BAM', 'SortedByCoordinate']) #BAM file
if run_pass == 2:
    other_params.append('--twopassMode Basic')
@jobs_limit(star_batch)
@follows(star_index)
@mkdir(fastqFiles,formatter(),'{path[0]}/f02_bam')
#@files(get_fq)
@transform(trim_reads,formatter('\.gz'),'f02_bam/{basename[0]}.bam')
@check_if_uptodate(check_file_exists)
def run_star(input_file,output_file):
    print input_file,output_file
    n = num_thread2use(star_batch,len(fastqFiles),thread)
    STAR(input_file,output_file,star_db,n,annotation,other_params)

#--------------------- 4. make tag_directory ------------------------------------------------------
@follows(run_star)
@mkdir(fastqFiles,formatter(),'{path[0]}/f03_tags')
#@check_if_uptodate(check_file_exists)
@transform(run_star,formatter('\.bam'),'f03_tags/{basename[0]}')
def make_tag(input_bam,out_dir):
    print(input_bam,out_dir)
    make_tag_directory(input_bam, out_dir, ref_fa)


#--------------------- 5. make bedgraph file for visualization ------------------------------------
@follows(make_tag)
@mkdir(fastqFiles,formatter(),'{path[0]}/bedgraph_files')
#@check_if_uptodate(check_file_exists)
@transform(make_tag,formatter('.'),['bedgraph_files/{basename[0]}_pos.bedgraph','bedgraph_files/{basename[0]}_neg.bedgraph'])
def make_bedgraph(in_dir,out_files):
    make_bedgraph_file(in_dir, out_files)


# #--------------------- 6. find peaks ------------------------------------------------------
def get_input_for_peak_call():
    gro_cap = [f for f in os.listdir('f03_tags') if 'mSTART' in f and 'input' not in f]
    gro_seq = [f for f in os.listdir('f03_tags') if 'input'  in f]
    comb = list(itertools.product(gro_cap,gro_seq))
    for com in comb:
        out = com[0] + '_and_' + com[1] + '_bg_' + str(peak_bg)
        yield ['f03_tags/'+f for f in com],'f04_peaks/' + out + '.peak'


@jobs_limit(thread)
@follows(make_bedgraph)
@mkdir(fastqFiles,formatter(),'{path[0]}/f04_peaks')
@files(get_input_for_peak_call)
#@check_if_uptodate(check_file_exists)
def find_peak(input_files,output_file):
    print input_files,output_file
    bg = '-F %d' % (peak_bg) #Enrichment over backgroun
    find_peaks(input_files[0],output_file,'tss',input_files[1],[bg])


#--------------------- 7. Merge peaks ------------------------------------------------------
@follows(find_peak)
@merge(find_peak,'f04_peaks/merge_bg_' + str(peak_bg) +'.peak')
#@check_if_uptodate(check_file_exists)
def merge_peak(input_files,output_file):
    merge_peaks(input_files,output_file,'given',type_merge)


#--------------------- 8. Create bed file ------------------------------------------------------
@follows(merge_peak)
@mkdir(fastqFiles,formatter(),'{path[0]}/bed_files')
@transform(merge_peak,formatter('\.peak'), 'bed_files/{basename[0]}')
#@check_if_uptodate(check_file_exists)
def create_bed_files(input_file,output_file):
    print input_file,output_file
    convert_peak_to_bed_file(input_file,output_file)


#--------------------- 9. Peak statistics ------------------------------------------------------
# @follows(create_bed_files)
# @mkdir(fastqFiles,formatter(),'{path[0]}/f05_statistics')
# @transform(merge_peak,formatter('\.peak'),'f05_statistics/{basename[0]}')
# #@check_if_uptodate(check_file_exists)
# def calculate_peak_statistics(input_file,output_file):
#     print input_file,output_file
#     wrap_stats_all_in_one(input_file,output_file,multiple_peaks_width,divergent_width,type_merge)

#--------------------- 10. annotate peaks ------------------------------------------------------
@jobs_limit(thread)
@follows(create_bed_files)
@mkdir(fastqFiles,formatter(),'{path[0]}/f06_annoPeaks')
@transform(merge_peak,formatter('\.peak'),'f06_annoPeaks/{basename[0]}.anno')
#@check_if_uptodate(check_file_exists)
def anno_peak(input_file,output_file):
    annotate_peaks(input_file,output_file,ref_fa,annotation)
    annotate_filter(output_file)
#--------------------- 11. histogram of peaks ------------------------------------------------------
@jobs_limit(thread)
@follows(anno_peak)
@mkdir(fastqFiles,formatter(),'{path[0]}/f07_histPeaks')
@transform(merge_peak,formatter('\.peak'),'f07_histPeaks/{basename[0]}.hist')
#@check_if_uptodate(check_file_exists)
def peak_cov_hist(input_file,output_file,pc=0): # input is peak file
    gro_cap = [f for f in os.listdir('f03_tags') if 'mSTART' in f and not 'input' in f]
    tag = ['f03_tags/' + t for t in gro_cap] # if t in input_file]
    hist(tag[0],output_file,ref_fa,annotation,mode='peak',peak=input_file,region=4000,res=25,pc=pc)
    hist_plot(output_file)


#--------------------- 12. histogram centered on annotation ------------------------------------------------------
@jobs_limit(thread)
@follows(peak_cov_hist)
@transform(merge_peak,formatter('\.peak'),'f07_histPeaks/{basename[0]}.hist_mrna')
#@check_if_uptodate(check_file_exists)
def peak_cov_hist_mrna(input_file,output_file,pc=0): # input is peak file
    gro_cap = [f for f in os.listdir('f03_tags') if 'mSTART' in f and not 'input' in f]
    tag = ['f03_tags/' + t for t in gro_cap] # if t in input_file]
    print('tag0',tag[0])
    print('input',input_file)
    print('output',output_file)
    #Peaks need to be filtered to peaks that have tags nearby (this helps the normalization and reduces computation)
    peakFileToPeakFile(mRNA_peak_file,input_file,distance=1000)
    hist(tag[0],output_file,ref_fa,annotation,mode='peak',peak=mRNA_peak_file+'filt',region=4000,res=25,pc=pc)
    hist_plot(output_file)
    heat_plot(output_file+'MatS',save_f = output_file + '_heat.png' )

#--------------------- 13. hist and descriptives separating the different annotations ------------------------------------------------------
@jobs_limit(thread)
@follows(peak_cov_hist_mrna)
@transform(anno_peak,formatter('.anno'),'f07_histPeaks/{basename[0]}_an.hist')
#@check_if_uptodate(check_file_exists)
def peak_cov_hist_anno(input_file,output_file,pc=0): # input is peak file
    print input_file,output_file
    gro_cap = [f for f in os.listdir('f03_tags') if 'mSTART' in f and not 'input' in f]
    tag = ['f03_tags/' + t for t in gro_cap]# if t in input_file]
    for mark in marks:
        hist(tag[0],output_file+'_' + mark,ref_fa,annotation,mode='peak',peak=input_file+ '_' + mark,region=4000,res=25,pc=pc)
        hist_plot(output_file+'_' + mark)



@follows(peak_cov_hist_anno)
def last_function():
    Message('START-seq analysis succeed',contact)
    
if __name__ == '__main__':
    try:
        #pipeline_printout(sys.stdout, [last_function], verbose=3)
        pipeline_run([last_function],forcedtorun_tasks = [last_function],multiprocess=thread,gnu_make_maximal_rebuild_mode = True,
                    touch_files_only=False,verbose=4)#,checksum_level=3)
    except:
        Message('START-seq analysis failed',contact)


'''
Generic pipeline for one file: 
STAR --genomeDir /data/genome/hamster/picr/picr_STAR_Db --readFilesCommand zcat --readFilesIn {input.fastq.gz} --runThreadN 6 --outFileNamePrefix f02_bam/{input}.bam --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --sjdbGTFfile /data/genome/hamster/picr/updated_final_sort.gff3 --sjdbGTFtagExonParentTranscript Parent
'''     