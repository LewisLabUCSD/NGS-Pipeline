import sarge,sys,os
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
mpl.style.use('ggplot')


def make_tag_directory(in_bam,tag_dir,ref_fa):
    '''make tag directory which extract mapping position into tsv file
    '''
    cmd = ('makeTagDirectory {o_dir} -genome {g} -checkGC \
            -single {bam}').format(o_dir=tag_dir,g=ref_fa,bam=in_bam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def rm_5GRO_ctrl(GRO5_tag,ctrl_tag):
    '''this function remove ctrl_tag coverage from GRO-Cap tag
    '''
    before_sub = GRO5_tag+'/genome.tages_before_sub.tsv'
    after_sub = GRO5_tag+'/genome.tags.tsv'
    os.rename(after_sub,before_sub)
    df_raw = pd.read_csv(before_sub,sep='\t',header=None)
    df_ctr = pd.read_csv(ctrl_tag+'/genome.tags.tsv',sep='\t',header=None)
    
    
def hist(tag_dir,hist_out,ref_fa,anno,mode='tss',peak='',region=4000,res=10,pc=3):
    '''this function gets tag coverage around tss
    * tag_dir: tag directory
    * anno: gff file or gtf file
    * pc: number of tags to consider at each position
    * region: length to conisder in the x axis. which means -2000 to 2000 around tss
    * res: resolution of the histogram.
    '''
    if anno.endswith('gtf'):
        anno = '-gtf ' + anno
    else:
        anno = '-gff ' + anno
    if mode == 'tss':
        cmd = ('annotatePeaks.pl tss {ref_fa} {anno} -size {size} -hist {bin} -d {dir} -pc {pc} > {out}').format(
                    ref_fa=ref_fa,anno=anno,size=str(region),bin=str(res),dir=tag_dir,pc=str(pc),out=hist_out)
    elif mode == 'peak':
        if peak == '':
            raise ValueError('input is empty')
        cmd = ('annotatePeaks.pl {peak} {ref_fa} {anno} -size {size} -hist {bin} -d {dir} -pc {pc} > {out}').format(
                    peak=peak,ref_fa=ref_fa,anno=anno,size=str(region),bin=str(res),dir=tag_dir,pc=str(pc),out=hist_out)
    sarge.run(cmd)
    
    
def hist_plot(hist_out):
    #Visualize histogram.
    plt.figure()
    df = pd.read_csv(hist_out,sep='\t',header=0,names=['Distance from TSS','Coverage','+ Tags','- Tags'])
    plt.plot(df['Distance from TSS'],df['+ Tags'],label='+ Tags')
    plt.plot(df['Distance from TSS'],df['- Tags'],label='- Tags')
    plt.xlim([-500,500])
    plt.xlabel('Distance from TSS')
    plt.ylabel('Reads per bp per TSS')
    plt.axvline(x=0,c='k')
    plt.legend(loc='upper right')
    
    plt.savefig(os.path.splitext(hist_out)[0]+'.png')
    
# if __name__ == '__main__':
#     import glob
#     hists = glob.glob('/data/shangzhong/TSS/fq/f03_tags/*/hist.txt')
#     for h in hists:
#         plt.figure()
#         hist_plot(h)


def find_peaks(tag_dir,out_file,peak_style,control_dir,otherParams=['']):
    '''find peaks
    '''
    cmd = 'findPeaks {tag} -style {style} -o {out} -i {control} {other}'.format(
                tag=tag_dir,style=peak_style,out=out_file,control=control_dir,
                other=' '.join(otherParams))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

    
def merge_peaks(input_files,output_file,dist):
    '''
    * input_files: a list of peak files, name format is 5gro_and_gro
    * otuput_file: final merged peak file
    '''
    cmd = ('mergePeaks -d {dist} {in_files} > {out}').format(dist=str(dist),
                        in_files=' '.join(input_files),out=output_file)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def annotate_peaks(peak_file,output_file,ref_fa,annotation):
    '''
    this function annotate peaks, basically get closes TSS to each peak.
    '''
    if annotation.endswith('gtf'):
        anno = '-gtf ' + annotation
    else:
        anno = '-gff ' + annotation
    cmd = 'annotatePeaks.pl {peaks} {genome} {annotation} > {out}'.format(peaks=peak_file,genome=ref_fa,
                            annotation=anno, out=output_file)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def filter_anno_peak(in_peak_file,filter_peak_file):
    '''this function extracts the reliable TSS from the peak file
    The rule is: for each 5GRO, get overlap of peaks against different GROseq.
    Then get union set from pevious peaks.
    '''
    df = pd.read_csv(in_peak_file,sep='\t',header=0)
    gro5 = []
    gro = []
    for f in df['Focus Ratio/Region Size']:
        files = f.split('|')
        for sub_f in files: # sub_f is peak file result
            peaks = sub_f.split('_and_')  # peaks has 5gro file and gro file name
            for p in peaks:
                if p in gro5 or p in gro:
                    continue
                else:
                    if '5GRO' in p:
                        gro5.append(p)
                    else:
                        gro.append(p)
    print gro5,gro
    def extract_peak(gro5,gro,fns):
        '''fns is the splited filename in 6th column of annotated peak file'''
        keep = []
        res = True
        for g5 in gro5:
            keep.append([g5+'_and_'+g in fns for g in gro])
        for k in keep:
            if False in k:
                res=False
        return res
        
    # filter out the annopeak
    cri = df['Focus Ratio/Region Size'].map(lambda x: extract_peak(gro5,gro,x))
    df = df[cri]
    df.to_csv(filter_peak_file,sep='\t',index=False)

# if __name__ == '__main__':
#     in_peak_file = '/data/shangzhong/TSS/fq/f05_annoPeaks/merge.anno'
#     filter_peak_file = '/data/shangzhong/TSS/fq/f05_annoPeaks/merge_filter.anno'
#     filter_anno_peak(in_peak_file,filter_peak_file)
                
            
    