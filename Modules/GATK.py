import sarge
import sys,re

def RealignerTargetCreator(dedupbam,interval,gatk,ref_fa,thread,gold_indels=['']):
    '''This function creates interval files for realigning.
    Input is deduplicated sorted bam files. reference is 
    fasta file.
    '''
    cmd = ('java -jar {gatk} -T RealignerTargetCreator '
           '-R {ref_fa} -I {dedup} -o {output} -nt {thread} ').format(
            gatk=gatk,ref_fa=ref_fa,dedup=dedupbam,output=interval,
            thread=str(thread))
    if gold_indels != ['']:
        gold_indels = ['-known ' + f for f in gold_indels]
        cmd = cmd + ' '.join(gold_indels)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def IndelRealigner(dedupBam,realiBam,gatk,ref_fa,interval,gold_indels=['']):
    '''This function realigns the deduped bam file to intervals
    reference is fasta file, target is target interval file.
    '''
    cmd = ('java -jar {gatk} -T IndelRealigner -R {ref_fa} '
           '-I {input} -targetIntervals {target} '
           '-o {output} ').format(gatk=gatk,ref_fa=ref_fa,
           input=dedupBam,target=interval,output=realiBam)
    if gold_indels != ['']:
        gold_indels = ['-known ' + f for f in gold_indels]
        cmd = cmd + gold_indels
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def HaplotypeCaller_DNA_gVCF(recalBam,vcf,gatk,ref_fa,thread,otherParameters=[]):
    '''
    this function does calling variant and stores the result 
    into the gVCF file.
    '''
    cmd = ('java -jar {gatk} -T HaplotypeCaller -R {ref_fa} -I {input} '
           '--emitRefConfidence GVCF -o {output} -nct {t}').format(
            gatk=gatk,ref_fa=ref_fa,input=recalBam,output=vcf,t=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def JointGenotype(raw_vcfs,gvcf,gatk,ref_fa,thread):
    ''' Merge all vcf files into one g.vcf file
    '''
    vcfs = ['--variant '+ f for f in raw_vcfs]
    cmd = ('java -Xmx100g -jar {gatk} -T GenotypeGVCFs -R {ref_fa} {vcf} '
           '-o {out} -nt {thread}').format(gatk=gatk,ref_fa=ref_fa,
            vcf=' '.join(vcfs),out=gvcf,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def SelectVariants(joint_variant,out_vcf,gatk,reference,extract_type,thread):
    """this function can extract either SNP or indel from the
    vcf file.
    """
    cmd = ('java -jar {gatk} -T SelectVariants -R {ref_fa} -V {input} '
           '-selectType {type} -o {output} -nt {thread}').format(gatk=gatk, 
            ref_fa=reference,input=joint_variant,type=extract_type,
            output=out_vcf,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def snpHardFilter(snp_vcf,out_vcf,gatk,ref_fa):
    """
    this function will filter the snps, output a gold standard snp database
    """
    filtercmd = ('QD < 3.0 || FS > 50.0 || MQ < 50.0 || HaplotypeScore > 10.0 '
                 '|| MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0')
    filtercmd = """'{filter}'""".format(filter=filtercmd)
    filtername = """'snp_filter'"""
    cmd = ('java -jar {gatk} -T VariantFiltration -R {ref_fa} -V {input} '
           '--filterExpression {filter} --filterName {filtername} '
           '-o {output}').format(gatk=gatk,ref_fa=ref_fa,input=snp_vcf,
                filter = filtercmd,filtername=filtername,output=out_vcf)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def indelHardFilter(indel_file,out_vcf,gatk,ref_fa):
    """
    this function filter the indels,output a gold standard indel database
    """
    filtercmd = ("QD < 2.0 || FS > 200.0 || ReadPosRankSum < -15.0")
    filtercmd = """'{filter}'""".format(filter=filtercmd)
    filtername = """'indel_filter'"""
    cmd = ('java -jar {gatk} -T VariantFiltration -R {ref_fa} -V {input} '
           '--filterExpression {filter} --filterName {filtername} '
           '-o {output}').format(gatk=gatk,ref_fa=ref_fa,input=indel_file,
            filter = filtercmd,filtername=filtername,output=out_vcf)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def HardFilter(raw_gvcf,gold_snp_indel,gatk,ref_fa,thread):
    """
    this function will apply artificial filter for snp and indel
    """
    snp_filter = re.sub('g\.vcf$','snp.vcf',raw_gvcf)
    indel_filter = re.sub('g\.vcf$','indel.vcf')
    SelectVariants(raw_gvcf,snp_filter,gatk,ref_fa,'SNP',str(thread))
    SelectVariants(raw_gvcf,indel_filter,gatk,ref_fa,'INDEL',str(thread))
    snpHardFilter(snp_filter,gold_snp_indel[0],gatk,ref_fa)
    indelHardFilter(indel_filter,gold_snp_indel[1],gatk,ref_fa)
    

def BaseRecalibrator_1(realiBam,table,gold_pair,gatk,ref_fa,thread):
    '''Step 1 of base recalibration.
    '''
    cmd = ('java -jar {gatk} -T BaseRecalibrator -R {ref_fa} '
           '-I {realignbam} -knownSites {snp} -knownSites {indel} '
           '-o {output} -nct {t}').format(gatk=gatk,ref_fa=ref_fa,
                realignbam=realiBam,snp=gold_pair[0],indel=gold_pair[1],
                output=table,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def BaseRecalibrator_2(realiBam,post_table,table,gold_pair,gatk,ref_fa,thread):
    '''Step 2 of base recalibration: get post table'''
    cmd = ('java -jar {gatk} -T BaseRecalibrator -R {ref_fa} '
           '-I {realignbam} -knownSites {snp} -knownSites {indel} -BQSR {table} '
           '-o {output} -nct {thread} && ').format(gatk=gatk,ref_fa=ref_fa,
            realignbam=realiBam,snp=gold_pair[0],indel=gold_pair[1],output=post_table,table=table,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
           

def BaseRecalibrator_3(table,plot,post_table,gatk,ref_fa):
    '''Step 3 of base recalibration: compare table and post table
    '''
    cmd = ('java -jar {gatk} -T AnalyzeCovariates -R {ref_fa} '
           '-before {table} -after {post_table} -plots {output}').format(
            gatk=gatk,ref_fa=ref_fa,table=table,post_table=post_table,output=plot)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def BaseRecalibrator_4(realiBam,recalBam,gatk,ref_fa,gold_pair,table,thread):
    ''' Step 4 of base recalibration: recalibrate the base quality.
    '''
    cmd = ('java -jar {gatk} -T PrintReads -R {ref_fa} -I {input} -BQSR {table} '
           '-o {output} -nct {thread}').format(gatk=gatk,
            ref_fa=ref_fa,input=realiBam,table=table,output=recalBam,
            thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def CombineSNPandINDEL(vcfFiles,outvcf,gatk,ref_fa,otherParams=[]):
    """
    This function combines the vcf files.
    * gatk: gatk software pathway
    * ref_fa: reference genome fasta file
    * variantFiles: a list of vcf files that need to be combined
    * argus: additional argument
    """
    variCmd = ' '.join(['-V '+vcf for vcf in vcfFiles])
    other = ' '.join(otherParams)
    cmd = ('java -jar {gatk} -R {ref_fa} -T CombineVariants '
           '{varis} -o {outputVcf} {other}').format(gatk=gatk,ref_fa=ref_fa,
            varis=variCmd,outputVcf=outvcf,argu=other)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


#===============================================================================
#                     RNA part functions
#===============================================================================
def splitN(dedupBam,splitBam,gatk,ref_fa):
    '''This function splits reads due to wrong splicng by STAR'''
    cmd = ('java -jar {gatk} -T SplitNCigarReads -R {ref_fa} '
            '-I {input} -o {output} -rf ReassignOneMappingQuality '
            '-RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS').format(
            gatk=gatk,ref_fa=ref_fa,input=dedupBam,output=splitBam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def HaplotypeCaller_RNA_VCF(recalBam,vcf,gatk,ref_fa,thread='1'):
    """
    This function calls variants in RNAseq
    """
    cmd = ('java -jar {gatk} -T HaplotypeCaller -R {ref_fa} '
    '-I {input} -dontUseSoftClippedBases ' 
    '-stand_call_conf 20.0 -stand_emit_conf 20.0 -o {output} -nct {thread}').format(
    gatk=gatk,ref_fa=ref_fa,input=recalBam,output=vcf,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def RNA_Vari_Filter(vcf,filterVCF,gatk,ref_fa):
    """
    This function filter out the results of the vari call 
    """
    FS = """'FS > 30.0'"""
    QD = """'QD < 2.0'"""
    cmd = ('java -jar {gatk} -T VariantFiltration -R {ref_fa} '
            '-V {input} -window 35 -cluster 3 -filterName FS '
            '-filter {FS} -filterName QD -filter {QD} '
            '-o {output}').format(gatk=gatk,ref_fa=ref_fa,
                input=vcf,FS=FS,QD=QD,output=filterVCF)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    

def RNA_BaseRecalibrator_1(realiBam,table,gatk,ref_fa,gold_vcf,thread):
    '''step 1 of base recalibration,generate a table'''
    cmd = ('java -jar {gatk} -T BaseRecalibrator -R {ref_fa} '
            '-I {realignbam} -knownSites {gold} '
            '-o {output} -nct {thread}').format(gatk=gatk,ref_fa=ref_fa,
            realignbam=realiBam,gold=gold_vcf,output=table,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def RNA_BaseRecalibrator_2(realiBam,post_table,table,gatk,ref_fa,gold_vcf,thread='1'):
    '''Step 2 of base recalibration,generate post table'''
    cmd = ('java -jar {gatk} -T BaseRecalibrator -R {ref_fa} '
           '-I {realignbam} -knownSites {gold} -BQSR {table} '
           '-o {output} -nct {thread}').format(gatk=gatk,
            ref_fa=ref_fa,realignbam=realiBam,gold=gold_vcf,
            output=post_table,table=table,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


def RNA_BaseRecalibrator3(table,plot,post_table,gatk,ref_fa):
    '''Step 3 of base recalibration, compare the two tables'''
    cmd = ('java -jar {gatk} -T AnalyzeCovariates -R {ref_fa} '
           '-before {table} -after {post_table} -plots {output}').format(
            gatk=gatk,ref_fa=ref_fa,table=table,post_table=post_table,output=plot)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    
    
def RNA_BaseRecalibrator4(realiBam,recalBam,gatk,table,ref_fa,gold_vcf,thread='1'):    
    '''Step 4 of base recalibration'''
    cmd = ('java -jar {gatk} -T PrintReads -R {ref_fa} '
           '-I {input} -BQSR {table} -o {output} -nct {thread}').format(gatk=gatk,
            ref_fa=ref_fa,input=realiBam,table=table,output=recalBam,thread=str(thread))
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)
    
    
