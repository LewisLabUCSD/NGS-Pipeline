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

if p.rna96_seq_df.endswith('csv'):
    rna96_seq_df = pd.read_csv(p.rna96_seq_df)


    # used_seq =  #'exp' # alternative: theor
    seq_postfix  = '_Johan96_' + p.use_seq
    # out_file = "test.gff"

    in_seq =[]
    for i, (Cultivation, human_symbol, seq) in rna96_seq_df[['Cultivation','human_symbol', 'Construct %s seq (nt)'%(p.use_seq) ]].iterrows():
        in_seq.append(SeqRecord(seq=Seq(seq), id = Cultivation, name = human_symbol + seq_postfix, description= 'from_johan_'))
elif p.rna96_seq_df.endswith('.fasta') | p.rna96_seq_df.endswith('.fa'):
    in_seq = list(SeqIO.parse(p.rna96_seq_df, "fasta"))
else:
    raise ValueError

start_end = {}
N_length = 0
end_ = 0
i_gene_start = 28978 + 1 # NC_007936.1     RefSeq  gene    14123   15265   .       +       .       ID=gene28978;Dbxref=GeneID:3979182;Name=CYTB;gbkey=Gene;gene=CYTB;gene_biotype=protein_coding
i_rna_start = 41557 + 1
i_id_start = 534773 + 1
recs = []
for i, seq in enumerate(in_seq):
    new_seq = ''
    # i+=i_start
    # new_seq+= 'N'*N_length ## append
    start_ = len(new_seq)#N_length + end_
    new_seq+= seq.seq
    end_ = len(new_seq) #start_ + len(seq)
    start_end[i] = (start_, end_)

    sub_features = []
    gene_feature = SeqFeature(FeatureLocation(start_, end_),  type = 'gene', strand=1,
                              qualifiers={"source": p.seq_source,
                                          "ID": "gene%i"%(i + i_gene_start),
                                          'gbkey': 'Gene',
                                          'gene': seq.name,
                                          'Name': seq.name})
    mrna_feature = SeqFeature(FeatureLocation(start_, end_),  type = 'mRNA', strand=1,
                              qualifiers={"source": p.seq_source,
                                          "ID": "rna%i"%(i + i_rna_start),
                                          'gbkey': 'mRNA',
                                          'gene': seq.name,
                                          'Name': seq.name + '_mRNA'}) ## all exons
    exon_feature = SeqFeature(FeatureLocation(start_, end_),  type = 'exon', strand=1,
                              qualifiers={"source": p.seq_source,
                                          "ID": "id%i"%(i + i_id_start), ## no name features for EXONs?
                                          'gbkey': 'mRNA',
                                          'gene': seq.name}) ## all exons
    mrna_feature.sub_features = [exon_feature]
    gene_feature.sub_features  = [mrna_feature]

    sub_features.append(gene_feature)
    top_feature = SeqFeature(FeatureLocation(0, end_), type="region", strand=1,
                             qualifiers={"source": p.seq_source,
                                         # "score": 10.0,
                                         # "other": ["Some", "annotations"],
                                         "ID": "RNA96_REF%i_0329"%(i)})
    top_feature.sub_features = sub_features
    rec = SeqRecord(new_seq, 'gene%i__'%(i) + seq.name , description=p.description)
    rec.features = [top_feature]
    recs.append(rec)


# seq = Seq("GATCGATCGATCGATCGATC")


with open('{name}_{use_seq}.gff'.format(name = p.new_scaffold_name,
                                        use_seq = p.use_seq ),
          "w") as out_handle:
    GFF.write(recs, out_handle)

SeqIO.write(recs,
            '{name}_{use_seq}.fa'.format(name = p.new_scaffold_name,
                                         use_seq = p.use_seq),
            'fasta')



cmd = 'cat {cho_assembly_fa} {generated96_fa} > {combined_fa}'.format(cho_assembly_fa = p.cho_assembly_fa,
                                                                     generated96_fa = '{name}_{use_seq}.fa'.format(name = p.new_scaffold_name,
                                                                                                                   use_seq = p.use_seq),
                                                                     combined_fa = p.ref_fa.format(use_seq = p.use_seq))
sarge.run(cmd)

cmd = 'cat {cho_assembly_gff} {generated96_gff} > {combined_gff}'.format(cho_assembly_gff = p.cho_assembly_gff,
                                                                         generated96_gff = '{name}_{use_seq}.gff'.format(name = p.new_scaffold_name,
                                                                                                                         use_seq = p.use_seq),
                                                                         combined_gff = p.gff.format(use_seq = p.use_seq))
sarge.run(cmd)
## cmd: cat chok1.fa ~/GoogleDrive/NGS-Pipeline/chr_RNA96.fa > chok1_96.fa
