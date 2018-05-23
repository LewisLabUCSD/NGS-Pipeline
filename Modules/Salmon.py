import sarge, re
import pandas as pd

def salmon_index(db_path,ref_fa):
    '''build index for salmon'''
    cmd = ('salmon index -t {rna} -i {idx}').format(rna=ref_fa,idx=db_path)
    print(cmd)
    sarge.run(cmd)


def salmon(input_files,out_path,salmon_index,thread,lib=''):
    '''run salmon'''
    if lib == '':
        library = '-l A'
    else:
        library = '-l ' + lib
    if len(input_files) == 2:
        cmd=('salmon quant -i {index} {lib} -p {t} -1 {f1} -2 {f2} -o {out_path}').format(
        index=salmon_index,lib=library,t=str(thread),f1=input_files[0],f2=input_files[1],out_path=out_path)
    else:
        cmd=('salmon quant -i {index} {lib} -p {t} -r {f1} -o {out_path}').format(
        index=salmon_index,lib=library,t=str(thread),f1=input_files[0],out_path=out_path)
    print(cmd)
    sarge.run(cmd)

def get_gene_expression(quant, gff_fn, anno_type='ncbi'):
    '''this function combines rna expression into gene expression'''
    # 
    out = quant[:-2]+'gene.txt'
    df = pd.read_csv(quant,sep='\t',header=0)
    df = df[['Name','TPM','NumReads']]
    def extract_id(anno, anno_type):
        anno = anno[8]
        if anno_type == 'ncbi':
            rna = re.search('(?<=transcript_id=).+?(?=$|[;,])',anno).group(0)
            gid = re.search('(?<=GeneID:).+?(?=$|[;,])',anno).group(0)
            gnm = re.search('(?<=gene=).+?(?=$|[;,])',anno).group(0)
        elif anno_type == 'ensembl':
            rna = re.search('(?<=transcript_id=).+?(?=$|[;,])',anno).group(0)
            gid = re.search('(?<=gene_id=).+?(?=$|[;,])',anno).group(0)
            gnm = re.search('(?<=gene_name=).+?(?=$|[;,])',anno).group(0)
        return pd.Series([rna, gid, gnm])
    
    # read gff file
    gff_df = pd.read_csv(gff_fn,sep='\t',header=None,comment='#')
    if anno_type == 'ncbi':
        gff_df = gff_df[gff_df[2].values == 'mRNA']
    elif anno_type == 'ensembl':
        gff_df = gff_df[gff_df[2].values == 'transcript']
    gff_df[['rnaid','gid','gnm']] = gff_df.apply(lambda x: extract_id(x,anno_type),axis=1)
    rna_gid_dic = gff_df.set_index('rnaid')['gid'].to_dict()
    rna_gnm_dic = gff_df.set_index('rnaid')['gnm'].to_dict()
    df['geneid'] = df['Name'].map(lambda x: rna_gid_dic[x] if x in rna_gid_dic else x)
    df['gename'] = df['Name'].map(lambda x: rna_gnm_dic[x] if x in rna_gid_dic else x)
    df = df.groupby(['geneid','gename']).sum()
    df.to_csv(out,sep='\t')



