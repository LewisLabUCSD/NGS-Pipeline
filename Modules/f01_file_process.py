import os
from natsort import natsorted
import sarge
import re
import pandas as pd



class dic2obj:
    def __init__(self, **entries): 
        self.__dict__.update(entries)

def remove(files):
    """
    this function can remove files provided
    Arguments:  1. files: a list of files to be removed
    
    files: a list of files to be removed. [f1,f2,f3,...] or [[f1,f2],[f3],...] or with any depth of list layers
    """
    if isinstance(files,str):
        os.remove(files)
        try:
            os.remove(files + '.bai')
        except:
            pass
    if isinstance(files,list):
        for f in files:
            remove(f)
            try:
                remove(f+'.bai')
            except:
                continue

def check_file_exists(input_file, output_file):
    if not os.path.exists(output_file):
        return True, "Missing file %s" % output_file
    else:
        return False, "File %s exists" % output_file

def list_fq_files(file_path):
    """
    This function list all fastq files into a list
    """
    allFiles = [f for f in os.listdir(file_path) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]
    allFiles = natsorted(allFiles)
    fastqFiles = []  # this list is going to stroe the paired or single file for running aligner
    while len(allFiles) > 1:           # this is to append the single end or pair end files into a list.
        if allFiles[0].endswith(".fastq.gz"):
            index = allFiles[0].index(".fastq.gz")
            if allFiles[1][index-2:index] == '_2':
                fastqFiles.append(allFiles[:2])
                del allFiles[:2]
            else:
                fastqFiles.append(allFiles[:1])
                del allFiles[:1]
        
        if len(allFiles) != 0:        
            if allFiles[0].endswith(".fq.gz"):
                index = allFiles[0].index(".fq.gz")
                if allFiles[1][index-2:index] == '_2':
                    fastqFiles.append(allFiles[:2])
                    del allFiles[:2]
                else:
                    fastqFiles.append(allFiles[:1])
                    del allFiles[:1]
    if len(allFiles) == 1:
        fastqFiles.append(allFiles)

    return fastqFiles


def replace_filename(inputfile,input_pattern,out_pattern):
    """
    This function generates the outputfile name to address a problem that transform failed to do:
       make the outputfile the same length with inputfile when inputfile length varies
    * inputfile: list. [[f1.fq.gz]...] or [[f1.fq.gz,f2.fq.gz]]
    * input_pattern: a pattern in input file.
    * out_pattern: a pattern for output file
    """
    outFile = []
    regex = re.compile(input_pattern)
    for infile in inputfile:
        res = []
        for fn in infile:
            output = regex.sub(out_pattern,fn)
            res.append(output)
        outFile.append(res)
    return inputfile,outFile
# result = replace_filename([['f_1.fq.gz','f_2.fq.gz']],'^','')
# print result

def num_thread2use(jobs_per_batch,len_of_jobs,given_thread):
    """
    This function calculates how many thread to use for each job given the jobs to run per batch and total number of jobs
    Some times total number of jobs is less than provided, in this case we can assign more thread to each job to run fast.
    """
    jobs = min(jobs_per_batch,len_of_jobs)
    if jobs ==0:
        thread = 1
    else:
        thread = int(given_thread/jobs)
    if thread == 0:
        thread = 1
    return thread


def Message(string,email):
    """
    This function send message to email when it run. 
    Used to calculate the time code runs.
    """
    cmd = ('echo {quote}|mailx -s "{string}" {email}').format(quote="",string=string,email=email)
    sarge.run(cmd)


def id_symbol_conversion(input_file,output_file,gene2refseq,tax_id,sym2ID='yes'):
    """This function convers count file based on gene symbol to gene id
    * inputfile: 2 columns. ['symbol','count']
    * outputfile: 2 columns. ['geneid','count']"""
    # 1. build {symbol:id conversion}
    df = pd.read_csv(gene2refseq,sep='\t',header=None,usecols=[0,1,15],names=['tax','geneid','symbol'],comment='#',compression='gzip')
    df = df[df['tax'].values==int(tax_id)]
    sym_id_dict = df.set_index('symbol')['geneid'].to_dict()
    # 2. transfer symbol -> id
    symbol_df = pd.read_csv(input_file,sep='\t',header=None,names=['symbol','count'])
    symbol_df['geneid'] = symbol_df['symbol'].map(lambda x: sym_id_dict[x] if x in sym_id_dict else x)
    # 3. output
    symbol_df[['geneid','count']].to_csv(output_file,sep='\t',header=None,index=False)
    os.remove(input_file)

     