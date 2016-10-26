from Modules.GeneMark import geneMark_ES
import os,sarge
from Bio import SeqIO
import glob
from natsort import natsorted
import multiprocessing as mp
import sys
# database files
ref_fa = '/data/genome/hamster/multi_pacbio_assemble/picr.fa'
rna_fa = '/data/shangzhong/Picr_assembly/Annotation/hamster_rna.fa'
refseq_pr = '/data/shangzhong/Picr_assembly/Annotation/hamster_pr.fa'
# pathways
path = '/data/shangzhong/Picr_assembly/Annotation'
genemark_path = path + '/genemark'
exonerate_path = path + '/exonerate'
PASA_path = path + '/PASA'
thread = '20'
# PASA parameters
pasa = '/home/shangzhong/Installation/PASApipeline-2.0.2'
ppl_fn = pasa + '/scripts/Launch_PASA_pipeline.pl'
config = pasa + '/pasa_conf/pasa.alignAssembly.Template.txt'
cmp_config = pasa + '/pasa_conf/pasa.annotationCompare.Template.txt'
load_fn = pasa + '/scripts/Load_Current_Gene_Annotations.dbi'
gff3_validate_fn = pasa + '/misc_utilities/pasa_gff3_validator.pl'
# evm weight file
evm_path = path + '/EVM'
exon2align_gff = '/home/shangzhong/Installation/EVidenceModeler-1.1.1/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl'

#===============================================================================
#                     1. run GeneMark
#===============================================================================
os.chdir(genemark_path)
# genemark_gff = geneMark_ES(ref_fa)
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
    for record in handle:
        SeqIO.write(record,open(path+'/file'+str(file_n)+'.fa','a'),'fasta')
        pr_n += 1
        if pr_n % int(item_per_file) == 0:
            file_n +=1

def exonerate2gff(gffs,out_gff):
    '''This function transfer exonerate gff file to standard gff format.
    gffs: a list of gff files
    out_gff: output final gff to store information
    '''
    out_handle = open(out_gff,'w')
    n = 1
    for gff in gffs:
        cds = []
        for line in open(gff):
            if line.startswith('#') or line.startswith('Command') or line.startswith('Hostname') or line.startswith(' ') or line.startswith('--'):
                continue
            else:
                item = line.strip().split('\t')
                if item[2] == 'cds':
                    cds.append(line.strip().split('\t'))
                elif item[2] == 'similarity':
                    info = item[8].split(';')
                    pr = info[1].split()[1]
                    length = 0
                    start = 1; end = 1
                    for c in cds:
                        length += int(c[4]) - int(c[3]) + 1
                        if length % 3 == 0:
                            end = length/3 
                            new_s = end + 1
                        else:
                            end = length/3 + 1
                            new_s = end
                        c[1] = 'exonerate'
                        c[2] = 'cds_match'
                        c.append(('ID=pr_{n};Target={pr} {s} {e}').format(n=n,pr=pr,s=start,e=end))
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
    # 1) split file
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
    
main_exonerate(ref_fa,refseq_pr,exonerate_path,thread,exon2align_gff,12)
    
#===============================================================================
#                     3. PASA Alignment assembly
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

# # 1. alignment assembly using gmap
# align_assemble(ppl_fn,config,ref_fa,rna_fa,thread) # 

# # 2. check gff compatability
# check_gff_compat(gff_fn,ppl_fn)
# # 3. load the gff file
# load_gff(gff_fn,ref_fa,load_fn,config)
# # 4. compare and update
# com_update(ref_fa,ppl_fn,cmp_config,rna_fa,thread)
#===============================================================================
#                     4. run EVM
#===============================================================================
def evm_partition(ref_fa,evm,gffs=[''],otherParams=['']):
    '''run evm to merge all the gff files'''
    cmd = ('{evm} --genome {ref} {gffs} {other} --segmentSize 50000000 \
     --overlapSize 10000 --partition_listing partitions_list.out').format(evm=evm,ref=ref_fa,
                                        gffs=' '.join(gffs),other=' '.join(otherParams))
    print(cmd)
    sarge.run(cmd)
    
def evm_cmd_list(out_fn,cmd_fn,evm,ref_fa,weight_fn,partition,gffs=['']):
    '''create cmd list for evm'''
    cmd = ('{evm} --genome {ref} --weights {w} {gffs} --output_file_name {out_fn} \
     --partitions {par} >  {cmd_l}').format(evm=evm,ref=ref_fa,
        w=weight_fn,gffs=' '.join(gffs),out_fn=out_fn,par=partition,cmd_l=cmd_fn)
    print(cmd)
    sarge.run(cmd)

def run_cmd(cmd):
    try:
        print(cmd);sys.stdout.flush()
        sarge.run(cmd)
    except:
        print cmd,'error'
        assert False

'''
evm = '/home/shangzhong/Installation/EVidenceModeler-1.1.1/EvmUtils'
os.chdir(evm_path)
genemark_gff = '/data/shangzhong/Picr_assembly/Annotation/genemark/genemark.gff.gff'
tr_gff = '/data/shangzhong/Picr_assembly/Annotation/PASA/picr_db.pasa_assemblies.gff3'#'/data/shangzhong/Picr_assembly/Annotation/PASA/picr_db.valid_gmap_alignments.gff3'
pr_gff = '/data/shangzhong/Picr_assembly/Annotation/exonerate/exonerate.gff'
evm_gffs = ['--gene_predictions '+genemark_gff,'--transcript_alignments '+tr_gff,'--protein_alignments '+pr_gff]
# 1. partition input
evm_partition(ref_fa,evm+'/partition_EVM_inputs.pl',evm_gffs)
# 2. generate command lines
evm_cmd_out = 'evm.out'
cmd_fn = 'commands.list'
weight_fn = '/home/shangzhong/Installation/EVidenceModeler-1.1.1/weights.txt'
evm_cmd_list(evm_cmd_out,cmd_fn,evm+'/write_EVM_commands.pl',ref_fa,weight_fn,'partitions_list.out',evm_gffs)
# 3. run commands
thread = 29
pool = mp.Pool(processes=int(thread))
cmds = open(cmd_fn).readlines()
for cmd in cmds:
    pool.apply_async(run_cmd,args=(cmd,))
pool.close()
pool.join()
'''










