import sarge
import re
import os

def Honey_pie(sortBam,sortTailBam,ref_fa,thread,tmp,otherParams=['']):
    """Honey pip extract soft clip reads and remap them
    """
    tailBam = re.sub('\.final\.bam$','.tail.bam',sortTailBam)
    cmd = ('Honey.py pie -o {tail} -n {thread} {input} {ref} --temp {tmp}').format(tail=tailBam,
                            thread=str(thread),input=sortBam,ref=ref_fa,tmp=tmp)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)
    # sort
    cmd = ('samtools sort -m 4G -@ {thread} -T {pre} -o {sortBam} {bam} ').format(
            thread=str(thread),pre=tailBam[:-4],sortBam=sortTailBam,bam=tailBam)
    print(cmd)
    sarge.run(cmd)
    # index
    cmd = ('samtools index {out} ').format(out=sortTailBam)
    print(cmd)
    sarge.run(cmd)
#     os.remove(sortBam)
    

def Honey_tails(finalBam,bamTail,otherParams=['']):
    """This function run Honey tail,culster the soft clipped reads
    """
    cmd = ('Honey.py tails -o {out} {input} ').format(input=finalBam,out=bamTail)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)


def Honey_spots(finalBam,spotFile,ref_fa,thread,otherParams=['']):
    """This function run Honey sorts.
    """
    cmd = ('Honey.py spots --reference {ref} -n {thread} -o {out} {input} ').format(
            input=finalBam,ref=ref_fa,thread=str(thread),out=spotFile)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)
