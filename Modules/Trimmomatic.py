import sarge
import os

def Trimmomatic(fqFiles,trim_fqFiles,trimmomatic,thread,phred='33',adapter_file=''):
    """This function run trimmomatic to trim reads"""
    # main parameters
    unpair = [f + 'unpair' for f in fqFiles]
    if len(fqFiles) == 1:
        trimCmd1st = ('java -jar {trim} SE -threads {thread} -phred{type} '
                              '{input} {output} ').format(trim=trimmomatic,thread = int(thread),
                            input = fqFiles,output=trim_fqFiles[0],type=phred)
        trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:5 TRAILING:3 MINLEN:22 TOPHRED33'
    elif len(fqFiles) == 2:
        trimCmd1st = ('java -jar {trim} PE -threads {thread} -phred{type} {fastq1} {fastq2} '
                '{Trimmed1} {unpair1} {Trimmed2} {unpair2} ').format(trim=trimmomatic,
                    thread=int(thread),type=phred,fastq1 = fqFiles[0], fastq2=fqFiles[1], 
                    Trimmed1 = trim_fqFiles[0], Trimmed2 = trim_fqFiles[1],unpair1=unpair[0],unpair2=unpair[1])
        trimCmd2nd = 'SLIDINGWINDOW:5:10 LEADING:15 TRAILING:10 MINLEN:36 TOPHRED33'
    # adapter file
    if adapter_file != '':
        adaptCmd = 'ILLUMINACLIP:{adapter}:2:30:10 '.format(adapter=adapter_file)
    else:
        adaptCmd = ''
    cmd = trimCmd1st + adaptCmd + trimCmd2nd
    print(cmd)
    sarge.run(cmd)
    os.remove(unpair[0]);os.remove(unpair[1])
