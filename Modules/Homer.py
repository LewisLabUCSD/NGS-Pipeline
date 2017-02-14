import sarge,sys


def make_tag_directory(in_bam,tag_dir,ref_fa):
    '''make tag directory which extract mapping position into tsv file
    '''
    cmd = ('makeTagDirecotry {o_dir} -genome {g} -checkGC \
            -single {bam}').format(o_dir=tag_dir,g=ref_fa,bam=in_bam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def find_peaks(tag_dir,out_dir,peak_style,control_dir):
    '''find peaks
    '''
    cmd = 'findPeaks {tag} -style {style} -o {out} {control}'.format(
                tag=tag_dir,style=peak_style,out=out_dir,control=control_dir)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    