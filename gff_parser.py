from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sarge
import  pandas as pd


import yaml
import socket
from  Modules import f01_file_process
parameter_file =  'Parameters/RNAseq_count_%s.yaml'%socket.gethostname()
with open(parameter_file,'r') as f:
    doc = yaml.load(f)
p = f01_file_process.dic2obj(**doc)


rna96_seq_df = pd.read_csv(p.rna96_seq_df)


# used_seq =  #'exp' # alternative: theor
new_scaffold_name = "chr_RNA96"
seq_postfix  = '_Johan96'
# out_file = "test.gff"

in_seq =[]
for i, (Cultivation, human_symbol, seq) in rna96_seq_df[['Cultivation','human_symbol', 'Construct %s seq (nt)'%(p.use_seq) ]].iterrows():
    in_seq.append(SeqRecord(seq=Seq(seq), id = Cultivation, name = human_symbol + seq_postfix, description= 'from_johan_'))
# in_seq = list(SeqIO.parse("../ppi/output/96Seq_Cultivation_nt.fasta", "fasta"))
start_end = {}
N_length = 10000
end_ = 0
sub_features = []
new_seq = ''
i_gene_start = 28978 + 1 # NC_007936.1     RefSeq  gene    14123   15265   .       +       .       ID=gene28978;Dbxref=GeneID:3979182;Name=CYTB;gbkey=Gene;gene=CYTB;gene_biotype=protein_coding
i_rna_start = 41557 + 1
i_id_start = 534773 + 1
for i, seq in enumerate(in_seq):
    # i+=i_start
    new_seq+= 'N'*N_length ## append
    start_ = len(new_seq)#N_length + end_
    new_seq+= seq.seq
    end_ = len(new_seq) #start_ + len(seq)
    start_end[i] = (start_, end_)

    gene_feature = SeqFeature(FeatureLocation(start_, end_),  type = 'gene', strand=1,
                              qualifiers={"source": "RNA96_REF",
                                          "ID": "gene%i"%(i + i_gene_start),
                                          'gbkey': 'Gene',
                                          'gene': seq.name,
                                          'Name': seq.name})
    mrna_feature = SeqFeature(FeatureLocation(start_, end_),  type = 'mRNA', strand=1,
                              qualifiers={"source": "RNA96_REF",
                                          "ID": "rna%i"%(i + i_rna_start),
                                          'gbkey': 'mRNA',
                                          'gene': seq.name,
                                          'Name': seq.name + '_mRNA'}) ## all exons
    exon_feature = SeqFeature(FeatureLocation(start_, end_),  type = 'exon', strand=1,
                              qualifiers={"source": "RNA96_REF",
                                          "ID": "id%i"%(i + i_id_start), ## no name features for EXONs?
                                          'gbkey': 'mRNA',
                                          'gene': seq.name}) ## all exons
    mrna_feature.sub_features = [exon_feature]
    gene_feature.sub_features  = [mrna_feature]

    sub_features.append(gene_feature)
top_feature = SeqFeature(FeatureLocation(0, end_), type="region", strand=1,
                         qualifiers={"source": "RNA96_REF",
                                     # "score": 10.0,
                                     # "other": ["Some", "annotations"],
                                     "ID": "RNA96_REF_0329"})
top_feature.sub_features = sub_features


# seq = Seq("GATCGATCGATCGATCGATC")
rec = SeqRecord(new_seq, new_scaffold_name, description = 'Johan96_scaffold')
rec.features = [top_feature]

with open('%s.gff'%(new_scaffold_name), "w") as out_handle:
    GFF.write([rec], out_handle)

SeqIO.write(rec, '%s.fa'%(new_scaffold_name), 'fasta')



cmd = 'cat {cho_assembly_fa} {generated96_fa} > {combined_fa}'.format(cho_assembly_fa = p.cho_assembly_fa,
                                                                     generated96_fa = '%s.fa'%(new_scaffold_name),
                                                                     combined_fa = p.ref_fa.format(use_seq = p.use_seq))
sarge.run(cmd)

cmd = 'cat {cho_assembly_gff} {generated96_gff} > {combined_gff}'.format(cho_assembly_gff = p.cho_assembly_gff,
                                                                         generated96_gff = '%s.gff'%(new_scaffold_name),
                                                                         combined_gff = p.gff.format(use_seq = p.use_seq))
sarge.run(cmd)
## cmd: cat chok1.fa ~/GoogleDrive/NGS-Pipeline/chr_RNA96.fa > chok1_96.fa