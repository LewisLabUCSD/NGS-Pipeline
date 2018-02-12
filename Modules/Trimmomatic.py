import sarge
import os
import gzip
import sys

def get_phred_score(fq):
    """This function get phred score of fastq.gz file
    """
    score_found = False
    with gzip.open(fq,'rb') as f:
        n = 0
        for line in f:
            n = n + 1
            line = line.rstrip()
            if n%4 == 0:  # only get quality line
                vals = [ord(c) for c in line]
                lmin = min(vals);lmax=max(vals)
                if lmin <= 50.:
                    return '33'
                    score_found=True
                    break
                if lmax >= 83.:
                    score_found=True
                    return '64'
                    break
    if score_found == False:
        raise 'could not find the phred score, need to manually set'
    
def Trimmomatic(fqFiles,trim_fqFiles,trimmomatic,thread,adapter_file='',min_len=36):
    """This function run trimmomatic to trim reads"""
    # main parameters
    unpair = [f + 'unpair' for f in fqFiles]
    phred = get_phred_score(fqFiles[0])
    if len(fqFiles) == 1:
        trimCmd1st = ('java -jar {trim} SE -threads {thread} -phred{type} '
                              '{input} {output} ').format(trim=trimmomatic,thread = int(thread),
                            input = fqFiles[0],output=trim_fqFiles[0],type=phred)
        trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:15 TRAILING:10 MINLEN:{len} TOPHRED33 '.format(len=min_len)
    elif len(fqFiles) == 2:
        trimCmd1st = ('java -jar {trim} PE -threads {thread} -phred{type} {fastq1} {fastq2} '
                '{Trimmed1} {unpair1} {Trimmed2} {unpair2} ').format(trim=trimmomatic,
                    thread=int(thread),type=phred,fastq1 = fqFiles[0], fastq2=fqFiles[1], 
                    Trimmed1 = trim_fqFiles[0], Trimmed2 = trim_fqFiles[1],unpair1=unpair[0],unpair2=unpair[1])
        trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:15 TRAILING:10 MINLEN:{len} TOPHRED33 '.format(len=str(min_len))
    # adapter file
    if adapter_file != '':
        adaptCmd = 'ILLUMINACLIP:{adapter}:2:30:10 '.format(adapter=adapter_file)
    else:
        adaptCmd = ''
    cmd = trimCmd1st + adaptCmd + trimCmd2nd
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    for un in unpair:
        if os.path.exists(un):
            os.remove(un)

def conda_Trimmomatic(fqFiles,trim_fqFiles,thread,adapter_file='',min_len=36):
    """This function run trimmomatic to trim reads"""
    # main parameters
    unpair = [f + 'unpair' for f in fqFiles]
    phred = get_phred_score(fqFiles[0])
    if len(fqFiles) == 1:
        trimCmd1st = ('trimmomatic SE -threads {thread} -phred{type} '
                              '{input} {output} ').format(thread = int(thread),
                            input = fqFiles[0],output=trim_fqFiles[0],type=phred)
        trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:15 TRAILING:10 MINLEN:{len} TOPHRED33 '.format(len=min_len)
    elif len(fqFiles) == 2:
        trimCmd1st = ('trimmomatic PE -threads {thread} -phred{type} {fastq1} {fastq2} '
                '{Trimmed1} {unpair1} {Trimmed2} {unpair2} ').format(
                    thread=int(thread),type=phred,fastq1 = fqFiles[0], fastq2=fqFiles[1], 
                    Trimmed1 = trim_fqFiles[0], Trimmed2 = trim_fqFiles[1],unpair1=unpair[0],unpair2=unpair[1])
        trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:15 TRAILING:10 MINLEN:{len} TOPHRED33 '.format(len=str(min_len))
    # adapter file
    if adapter_file != '':
        adaptCmd = 'ILLUMINACLIP:{adapter}:2:30:10 '.format(adapter=adapter_file)
    else:
        adaptCmd = ''
    cmd = trimCmd1st + adaptCmd + trimCmd2nd
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    for un in unpair:
        if os.path.exists(un):
            os.remove(un)