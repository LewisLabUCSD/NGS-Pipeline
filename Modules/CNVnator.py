import sarge,sys



def cnv_extract_bam(in_bam,out_root,others=['']):
    '''
    extract read mapping from bam files
    '''
    cmd = ('cnvnator -root {out} -unique -tree {bam} {other}').format(out=out_root,bam=in_bam,
                                                        other=' '.join(others))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def cnv_generate_hist(in_root,chr_path,bin_win,others=['']):
    '''
    generating histogram
    '''
    # 2. get histogram
    cmd = ('cnvnator -root {root} -his {bin} -d {dir} {other}').format(root=in_root,
                                bin=str(bin_win),dir=chr_path,other=' '.join(others))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def cnv_statistics(in_root,bin_win,others=['']):
    cmd = ('cnvnator -root {root} -stat {bin} {other}').format(root=in_root,bin=str(bin_win),other=' '.join(others))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def cnv_partitioning(in_root,bin_win,others=['']):
    cmd = ('cnvnator -root {root} -partition {bin} {other}').format(root=in_root,bin=str(bin_win),other=' '.join(others))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    
    
def cnv_call(in_root,out,bin_win,others=['']):
    cmd = ('cnvnator -root {root} -call {bin} {other} > {out}').format(root=in_root,bin=str(bin_win),out=out,other=' '.join(others))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    
# output_file = '/data/shangzhong/Pacbio/CHOS_illu_DNA/cnv/cnv/merge.txt'
# chr_path = '/data/shangzhong/Pacbio/CHOS_illu_DNA/cnv/cnv/scaffold'
# bin_win = 100
# others = ['-chrom NW_006887432.1']
# root = output_file[:-3] + 'root'
# cnv_generate_hist(root,chr_path,bin_win,others)
# # 3
# cnv_statistics(root,bin_win,others)
# # 4
# cnv_partitioning(root,bin_win,others)
# # 5
# cnv_call(root,output_file,bin_win,others)
    

