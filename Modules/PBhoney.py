import sarge
import re


def Honey_pie(bamFile,finalBam,ref_fa,thread,otherParams=['']):
    """Honey pip extract soft clip reads and remap them
    """
    tailBam = re.sub('final\.bam$','tail.bam',finalBam)
    cmd = ('Honey.py pie -o {tail} -n {thread} {input} {ref}').format(tail=tailBam,
                            thread=str(thread),input=bamFile,ref=ref_fa)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)
    cmd = ('samtools merge {out} {bam} {tail} ').format(out=finalBam,
            bam=bamFile,tail=tailBam)
    print(cmd)
    sarge.run(cmd)
    


def Honey_tails(finalBam,bamTail,otherParams=['']):
    """This function run Honey tail,culster the soft clipped reads
    """
    cmd = ('Honey.py tails {input} -o {out} ').format(input=finalBam,out=bamTail)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)


def Honey_spots(finalBam,spotFile,thread,otherParams=['']):
    """This function run Honey sorts.
    """
    cmd = ('Honey.py spots {input} -n {thread} -o {out} ').format(
            input=finalBam,thread=str(thread),out=spotFile)
    cmd = cmd + ' '.join(otherParams)
    print(cmd)
    sarge.run(cmd)
