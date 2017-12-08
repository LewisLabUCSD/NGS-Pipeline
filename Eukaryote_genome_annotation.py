import os,re
from Bio import SeqIO
import sarge
import glob
from natsort import natsorted
import multiprocessing as mp
import pandas as pd
from Bio import Entrez
import sys
from Bio import Entrez
Entrez.email = 'shl198@eng.ucsd.edu'

# database files
ref_fa = '/data/genome/hamster/multi_pacbio_assemble/picr.fa'
rna_fa = '/data/shangzhong/Picr_assembly/Annotation/hamster_rna.fa'
refseq_pr = '/data/shangzhong/Picr_assembly/Annotation/hamster_pr.fa'
hamster_id = '/data/shangzhong/Database/hamster/hamster_all_id.txt'
# pathways
path = '/data/shangzhong/Picr_assembly/Annotation'
organism = 'hamster'
# exonerate parameters
exonerate_path = path + '/exonerate'
pr_gff = exonerate_path + '/exonerate.gff'
# PASA parameters
PASA_path = path + '/PASA'
pasa = '/home/shangzhong/Installation/PASApipeline-2.0.2'
ppl_fn = pasa + '/scripts/Launch_PASA_pipeline.pl'
config = pasa + '/pasa_conf/pasa.alignAssembly.Template.txt'
cmp_config = pasa + '/pasa_conf/pasa.annotationCompare.Template.txt'
load_fn = pasa + '/scripts/Load_Current_Gene_Annotations.dbi'
gff3_validate_fn = pasa + '/misc_utilities/pasa_gff3_validator.pl'
tr_gff = PASA_path + '/picr_db.pasa_assemblies.gff3'
#===============================================================================
#                     1. PASA Alignment assembly
#===============================================================================
def align_assemble(ppl_fn,config,ref_fa,rna_fa,thread,otherParameters=['']):
    '''This function do alignment assembly
    generate 4 type of files: 
    sample_mydb_pasa.assemblies.fasta :the PASA assemblies in FASTA format.
    sample_mydb_pasa.pasa_assemblies.gff3,.gtf,.bed :the PASA assembly structures.
    sample_mydb_pasa.pasa_alignment_assembly_building.ascii_illustrations.out :descriptions 
        of alignment assemblies and how they were constructed from the underlying transcript alignments.
    sample_mydb_pasa.pasa_assemblies_described.txt :tab-delimited format describing the contents
         of the PASA assemblies, including the identity of those transcripts that were assembled into the corresponding structure.
    '''
    cmd = ('{ppl} -c {config} -C -r -R -g {ref_fa} \
             -t {rna_fa} --ALIGNERS gmap --CPU {thread} {other}').format(ppl=ppl_fn,config=config,
                        ref_fa = ref_fa,rna_fa=rna_fa,thread=str(thread),other=' '.join(otherParameters))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)

def check_gff_compat(gff,ppl_fn,config):
    '''check the gff compatibility with pasa'''
    cmd = ('{ppl_fn} {gff}').format(ppl_fn=ppl_fn,gff=gff)
    sarge.run(cmd)

def load_gff(gff,ref_fa,ppl_fn,config):
    cmd = ('{ppl} -c {config} -g {ref} -P {gff}').format(ppl=ppl_fn,config=config,ref=ref_fa,gff=gff)
    print(cmd)
    sarge.run(cmd)

def com_update(ref_fa,ppl_fn,config,rna_fa,thread):
    '''compare the reads and update the annotation'''
    cmd = ('{ppl_fn} -c {config} -A -g {ref_fa} -t {rna} --CPU {t}').format(ppl_fn=ppl_fn,
                                        config=config,ref_fa=ref_fa,rna=rna_fa,t=str(thread))
    print(cmd)
    sarge.run(cmd)

def main_PASA(gff_fn,ppl_fn,config,ref_fa,rna_fa,thread):
    # 1. alignment assembly using gmap
    align_assemble(ppl_fn,config,ref_fa,rna_fa,thread) # 
    # 2. check gff compatability
    check_gff_compat(gff_fn,ppl_fn)
    # 3. load the gff file
    load_gff(gff_fn,ref_fa,load_fn,config)
    # 4. compare and update
    com_update(ref_fa,ppl_fn,cmp_config,rna_fa,thread)
#===============================================================================
#                     2. run exonerate
#===============================================================================
def exonerate(ref_fa,pr_fn,out_fn):
    '''map protein sequence to dna seq'''
    cmd = ('exonerate -m p2g -q {pr} -t {ref} --showalignment no \
    --showvulgar no --showtargetgff yes --minintron 20 --percent 50 \
    --score 100 --geneseed 250 -n 10 > {gff}').format(pr=pr_fn,ref=ref_fa,gff=out_fn)
    print(cmd)
    sarge.run(cmd)

def split_fa(fa,item_per_file,path):
    if not os.path.exists(path): os.mkdir(path)
    handle = SeqIO.parse(open(fa,'r'),'fasta')
    file_n = 0
    pr_n = 0
    out_fn = path+'/file'+str(file_n)+'.fa'
    if os.path.exists(out_fn): os.remove(out_fn)
    for record in handle:
        SeqIO.write(record,open(out_fn,'a'),'fasta')
        pr_n += 1
        if pr_n % int(item_per_file) == 0:
            file_n +=1
            out_fn = path+'/file'+str(file_n)+'.fa'
            if os.path.exists(out_fn): os.remove(out_fn)

def exonerate2gff(gffs,out_gff,g_type='evm'):
    '''This function transfer exonerate gff file to standard gff format.
    gffs: a list of gff files
    out_gff: output final gff to store information
    '''
    out_handle = open(out_gff,'w')
    n = 1
    m = 0
    for gff in gffs:
        cds = []
        for line in open(gff):
            if line.startswith('#') or line.startswith('Command') or line.startswith('Hostname') or line.startswith(' ') or line.startswith('--'):
                continue
            else:
                item = line.strip().split('\t')
                if item[2] == 'cds':
                    cds.append(line.strip().split('\t'))
                elif item[2] == 'gene' and g_type=='augustus':
                    item[1] = 'exonerate'
                    pr = item[8].split(';')[1].split(' ')[2]
                    item[8] = ('ID=gene_{n};Target={pr}').format(n=n,pr=pr)
                    out_handle.write('\t'.join(item) + '\n')
                elif item[2] == 'similarity':
                    info = item[8].split(';')
                    pr = info[1].split()[1]
                    length = 0
                    start = 1; end = 1
                    for c in cds:   # decide start of the AA of each exon
                        length += int(c[4]) - int(c[3]) + 1
                        if length % 3 == 0:
                            end = length/3 
                            new_s = end + 1
                        else:
                            end = length/3 + 1
                            new_s = end
                        c[1] = 'exonerate'
                        c[2] = 'cds_match'
                        m += 1
                        if g_type == 'evm':
                            m = n
                        c.append(('ID=pr_{m};Parent=gene_{n};Target={pr} {s} {e}').format(m=m,n=n,pr=pr,s=start,e=end))
                        start = new_s
                        out_handle.write('\t'.join(c) + '\n')
                    cds = []
                    n += 1
    out_handle.close()


def main_exonerate(ref_fa,refseq_pr,exonerate_path,thread,exon2align_gff,index_s=0,index_e=0):
    '''
    * refseq_pr: all protein seqeunces of the organism
    * path: path to store splited protein sequences.
    '''
    if not os.path.exists(exonerate_path): os.mkdir(exonerate_path)
    # 1) split protein fa file into many sub file, this is to parallel the process
    os.chdir(exonerate_path)
    if os.listdir(path) != []:
        split_fa(refseq_pr,100,exonerate_path)
    # 2) run exonerate for each file
    faFiles = natsorted(glob.glob('file*.fa'))
    if index_e == 0:
        faFiles = faFiles[index_s:]
    else:
        faFiles = faFiles[index_s:index_e]
    pool = mp.Pool(processes=int(thread))
    for f in faFiles:
        out = f[:-2]+'gff'
        pool.apply_async(exonerate,args=(ref_fa,f,out))
    pool.close()
    pool.join()
    # 3) merge the gff files
    exonerate_gff = 'exonerate.gff'
    if not os.path.exists(exonerate_gff):
        gff_fns = natsorted(glob.glob('file*.gff'))
        exonerate2gff(gff_fns,exonerate_gff)
    
# main_exonerate(ref_fa,refseq_pr,exonerate_path,thread,exon2align_gff)


#===============================================================================
#                     process the gmap results and exonerates results directly
#===============================================================================
#=============== 1. get all mapped geneid, rna_accession, pr_accession
def gene_rna_pr_id(hamster_id,gmap_gff,out_fn):
    '''this fnction get all gene rna pr id, including both refseq and gff information.
    * hamster_id: a file that has all ids in hamster.gff file
    * gmap_gff: gff results mapped using gmap
    * out_fn:  
    '''
    # rna accession in gff file
    ham_id_df = pd.read_csv(hamster_id,sep='\t',header=0)
    ham_id_df = ham_id_df.astype('str')
    ham_id_df['TrAccess'] = ham_id_df['TrAccess'].map(lambda x: x.split('.')[0])
    ham_id_df['PrAccess'] = ham_id_df['PrAccess'].map(lambda x: x.split('.')[0])
    rna_gene_dic = ham_id_df.set_index('TrAccess')['GeneID'].to_dict()
    rna_pr_dic = ham_id_df.set_index('TrAccess')['PrAccess'].to_dict()
    #-------- read rna gff file
    rna_df = pd.read_csv(gmap_gff,sep='\t',header=None,comment='#')
    # add rna accession column
    rna_df['rna_ac'] = rna_df[8].map(lambda x: re.search('(?<=ID=).+?(?=\.)',x).group(0))
    mrna = list(set(rna_df['rna_ac'].tolist()))
    # new rna in refseq compared to gff
    new_ref_rna = list(set(mrna) - set(rna_gene_dic.keys()))
    # get geneid for new ref_rna gene id
    for r in new_ref_rna:
        handle = Entrez.efetch(db='nucleotide',id=r,rettype='gb',retmode='text').read()
        geneid = re.search('(?<=GeneID:).+?(?=\")',handle).group(0)
        try:
            p = re.search('(?<=protein_id=\").+?(?=\.)',handle).group(0)
        except:
            p = '-'
        rna_gene_dic[r] = geneid
        rna_pr_dic[r] = p
    # transfer dic to dataframe
    r_g_df = pd.DataFrame.from_dict(rna_gene_dic,'index')
    r_g_df.columns = ['geneid']
    r_p_df = pd.DataFrame.from_dict(rna_pr_dic,'index')
    r_p_df.columns = ['pr_ac']
    g_r_p_df = pd.concat([r_g_df,r_p_df],axis=1)
    g_r_p_df['rna_ac'] = g_r_p_df.index
    g_r_p_df[['geneid','rna_ac','pr_ac']].to_csv(out_fn,sep='\t',index=False)

# gmap_exon_path = path + '/gmap_exonerate'
# if not os.path.exists(gmap_exon_path): os.mkdir(gmap_exon_path)
# os.chdir(gmap_exon_path)
# gmap_gff = PASA_path + '/gmap.spliced_alignments.gff3'
# g_r_p_id_fn = gmap_exon_path + '/01_gene_rna_pr.txt'
# gene_rna_pr_id(hamster_id,gmap_gff,g_r_p_id_fn)


def get_consensus_map(rna_df,pr_df,gene,rna_ac,pr_ac):
    '''this function check if the rna map and pr map have the same splice sites
    * rna_df: mRNA map to genome gff dataframe with additional rna_ac column
    * pr_df: protein map to genome dataframe with additional 'pr_ac' and 'pr_id' column
    '''
    if not rna_df.empty:
        # get rna scaffold name, if more than 1 scaffold then don't add it's annotation
        rna_chr = list(set(rna_df[0].tolist()))
        if len(rna_chr) != 1:
            assert False, rna_ac + ' map to multiple scaffolds'
        else:
            rna_chr = rna_chr[0]
        # get strand, if map to both strand don't output
        rna_str = list(set(rna_df[6].tolist()))
        if len(rna_str) != 1:
            assert False, rna_ac + ' map to both strands'
        else:
            rna_str = rna_str[0]
        # get rna splice sites
        rna_splice = natsorted(rna_df[3].tolist() + rna_df[4].tolist())
        # change exon id
        n = 1
        for i,row in rna_df.iterrows():
            item = row[8].split(';')
            iid = '.'.join(item[0].split('.')[:-1])
            anno = iid+' '+str(n)+';'+re.sub('Name.+?;','',';'.join(item[1:]))+';Parent='+rna_ac+';gene_id='+gene+';transcript_id='+rna_ac
#             anno = iid+'_'+str(n)+';'+ re.sub('Name','transcript_id',';'.join(item[1:]))+';Parent='+rna_ac+';gene_id='+gene
            rna_df.loc[i,8] = anno
            rna_df.loc[i,2] = 'exon'
            n += 1
    #--------------- process protein gff information
    if not pr_df.empty:
        pr_id = pr_df['pr_id'].tolist()[0]
        sub_pr_df = pr_df[(pr_df['pr_id'].values==pr_id) & (pr_df[0].values==rna_chr)].copy()
        # change cds id
        m = 1
        for i,row in sub_pr_df.iterrows():
            item = row[8].split(';')
            anno = 'ID='+pr_ac+'_'+str(m)+';'+';'.join(item[2:])+';protein_id='+pr_ac+';Parent='+rna_ac+';gene_id='+gene
            sub_pr_df.loc[i,8] = anno
            sub_pr_df.loc[i,2] = 'CDS'
            m += 1
        pr_splice = natsorted(sub_pr_df[3].tolist() + sub_pr_df[4].tolist())
        if sub_pr_df.shape[0] == 1:
            if not rna_splice[0]<pr_splice[0]<pr_splice[1]<rna_splice[1]:
                sub_pr_df = pd.DataFrame()
        else:
            rna_pr_sites_match = set(pr_splice[1:-1]).intersection(rna_splice)
            m_len = len(rna_pr_sites_match)
            pr_len = len(pr_splice[1:-1])
            if  m_len != pr_len:
                print pr_ac,m_len,'/',pr_len
            if len(pr_splice) > len(rna_splice):
                print 'protein has more splice than rna, rna/pr:',len(rna_splice),'/',len(pr_splice)
                sub_pr_df = pd.DataFrame()
    else:
        sub_pr_df = pr_df
    return rna_df,sub_pr_df,rna_chr,rna_splice[0],rna_splice[-1],rna_str


# import time
# process_start = time.time()

def gmap_exonerate_merge_gff(gmap_gff,exonerate_gff,gmap_exon_path,all_id_fn):
    #-------- read gmap gff file
    rna_df = pd.read_csv(gmap_gff,sep='\t',header=None,comment='#')
    rna_df['rna_ac'] = rna_df[8].map(lambda x: re.search('(?<=ID=).+?(?=\.)',x).group(0))
    # get multi mapping mRNAs
    multi_map_rna = list(set(rna_df[rna_df[8].map(lambda x: 'path2' in x)]['rna_ac'].tolist()))  
    # build gene rna protein id dictionary
    g_r_p_dic = {}
    g_r_p_id_fn = gmap_exon_path + '/01_gene_rna_pr.txt'
    handle = open(g_r_p_id_fn)
    for line in handle:
        item = line.strip().split('\t')
        if item[1] in multi_map_rna:
            continue
        if item[0] in g_r_p_dic:
            g_r_p_dic[item[0]][item[1]] = item[2]
        else:
            g_r_p_dic[item[0]] = {item[1]:item[2]}
    #-------- read exonerate gff file
    pr_df = pd.read_csv(pr_gff,sep='\t',header=None)
    pr_df['pr_ac'] = pr_df[8].map(lambda x: re.search('(?<=Target=).+?(?=\.)',x).group(0))
     
    def output_consensus_rna_pr(g,out_handle):
        '''this function finds the consistend rna and protein and out put to file
        g_r_p_dic: dictionary that has all the gene, rna and protein ids.
        g: gene id
        rna_df: rna gff dataframe
        pr_df: protein gff dataframe
        '''
        g_n = 0
        rna_pr_dic = g_r_p_dic[g]
        for rna in rna_pr_dic:
            pr = rna_pr_dic[rna]
            single_rna_df = rna_df[rna_df['rna_ac'].values==rna].copy()
            single_rna_df = single_rna_df.reset_index(drop=True)
            if not single_rna_df.empty:
                single_pr_df = pr_df[pr_df['pr_ac'].values==pr].copy()
                single_pr_df = single_pr_df.reset_index(drop=True)
                single_pr_df.loc[:,'pr_id'] = single_pr_df[8].map(lambda x: re.search('(?<=ID=).+?(?=;)',x).group(0))
                res_rna_df,res_pr_df,chrome,start,end,strand=get_consensus_map(single_rna_df,single_pr_df,str(g),rna,pr)
                if g_n == 0:
                    out_handle.write('\t'.join([chrome,'gmap_exonerate','gene',str(start),str(end),'.',\
                                                strand,'.','ID='+str(g)+';gene_id='+str(g)])+'\n')
                    g_n += 1
                if not res_rna_df.empty:
                    if rna.startswith('XR'):
                        feature = 'lncRNA'
                    else:
                        feature = 'mRNA'
                    out_handle.write('\t'.join([chrome,'gmap_exonerate',feature,str(start),str(end),'.',\
                                                strand,'.','ID='+rna+';Parent='+str(g)+';gene_id='+str(g)+';transcript_id='+rna])+'\n')
                    res_rna_df[range(9)].to_csv(out_handle,sep='\t',index=False,header=None)
                    if not res_pr_df.empty:
                        res_pr_df[range(9)].to_csv(out_handle,sep='\t',index=False,header=None)
     
    out_fn = '02_gmap_exonerate.gff'
    if os.path.exists(out_fn): os.remove(out_fn)
    with open(out_fn,'a') as f:
        for g in g_r_p_dic.keys():
            output_consensus_rna_pr(g,f)
    # define a function to find the start and end position of each gene
    def get_gene_s_e(gene_df):
        pos = gene_df[3].tolist() + gene_df[4].tolist()
        gene_df.iloc[0,3] = min(pos)
        gene_df.iloc[0,4] = max(pos)
        return gene_df
    # correct gene coordinates
    gff_df = pd.read_csv(out_fn,sep='\t',header=None)
    gff_df['geneid'] = gff_df[8].map(lambda x: re.search('(?<=gene_id=).+?(?=;|$)',x).group(0))
    res_df = gff_df.groupby('geneid').apply(get_gene_s_e)
    # add gene name 
    all_id_df = pd.read_csv(all_id_fn,sep='\t',header=0)
    all_id_df = all_id_df.astype('str')
    g_s_dic = all_id_df.set_index('GeneID')['GeneSymbol'].to_dict()
    res_df[8] = res_df.apply(lambda row: row[8]+';gene_name='+g_s_dic[row['geneid']] if row['geneid'] in g_s_dic else row[8]+';gene_name='+row['geneid'],axis=1)
    res_df[range(9)].to_csv('02_gmap_exonerate.gff',sep='\t',index=False,header=None)
    
# gmap_gff = PASA_path+'/gmap.spliced_alignments.gff3'
# gmap_exonerate_merge_gff(gmap_gff,pr_gff,gmap_exon_path,hamster_id)
# 
# print time.time() - process_start







