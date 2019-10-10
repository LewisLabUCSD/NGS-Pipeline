#!/bin/bash

# RNAseq qc and quantification
#QC: Picard and RNA-seQC
#quantification: HT-seq and RSME

##### Constants


##### Functions


usage()
{
    echo "usage: STARpipeline.sh -f1 fastq1 -f2 fastq2 -sid sample_id -starIndx star_index_directory -BAM STAR_bam_path -refFASTA refereceGenome -GFF referenceGTF -RSEM RSEMreference -ncores num_cores"
}


##### Main
#fastq1="/data/vahid/RNAseqB5STAR/1870-01_W4G6DG07_S55/fastq1.gz"
#fastq2="/data/vahid/RNAseqB5STAR/1870-01_W4G6DG07_S55/fastq2.gz"
#sample_id="1870-01_W4G6DG07_S55"
#bam_path="/data/vahid/RNAseqB5STAR/1870-01_W4G6DG07_S55/"
#refGenome="/data/vahid/RNAseqPipeLine/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
#refGTF="/data/vahid/RNAseqPipeLine/gencode.v26.GRCh38.genes.gtf"
#RSEMpath="/data/vahid/RNAseqPipeLine/rsemRef/rsem_reference"
    
ncores=10
while [ "$1" != "" ]; do
    case $1 in
        -f1 | --fastq1 )        shift
                                fastq1=$1
                                ;;
        -f2 | --fastq2 )        shift
                                fastq2=$1
                                ;;
        -sid | --sample_id )    shift
                                sample_id=$1
                                ;;
        -starIndx | --Star_index_directory )      shift
                                star_index=$1
                                ;;
        -BAM | --BAMpath )      shift
                                bam_path=$1
                                ;;
        -refFASTA | --refGenome )        shift
                                refGenome=$1
                                ;;
        -GFF | --refGFF )       shift
                                refGTF=$1
                                ;;
        -ncores | --num_cores )       shift
                                ncores=$1
                                ;;
        -RSEM | --RSEMreference )       shift
                                RSEMpath=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
    esac
    shift
done

#mkdir -p -- ${bam_path}

#cp -f "$fastq1" ${bam_path}fastq1.gz
#cp -f "$fastq2" ${bam_path}fastq2.gz

bam_file_name="Aligned.sortedByCoord.out.bam"

STAR --runMode alignReads --runThreadN ${ncores} --genomeDir ${star_index} \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--readFilesIn $fastq1 $fastq2 \
--readFilesCommand zcat \
--outFileNamePrefix ${bam_path} \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType WithinBAM SoftClip \
--chimMainSegmentMultNmax 1 \
--outSAMattributes NH HI AS nM NM ch --outSAMattrRGline ID:rg1 SM:sm1

#Marking the duplicates using picard
mkdir -p -- ${bam_path}markedDup
java -jar /data/vahid/RNAseqPipeLine/picard.jar \
        MarkDuplicates I=${bam_path}${bam_file_name} \
        O=${bam_path}markedDup/${bam_file_name} \
        M=${bam_path}/markedDup/sample_id.marked_dup_metrics.txt \
        ASSUME_SORT_ORDER=coordinate

#Running Picard QC
mkdir -p -- ${bam_path}/PicardQC
java -jar /data/vahid/RNAseqPipeLine/picard.jar CollectAlignmentSummaryMetrics \
R=${refGenome} \
INPUT=${bam_path}${bam_file_name} \
OUTPUT=${bam_path}PicardQC/picard_QC.txt

#Running Picard fragment size
java -jar /data/vahid/RNAseqPipeLine/picard.jar CollectInsertSizeMetrics \
      I=${bam_path}${bam_file_name} \
      O=${bam_path}PicardQC/insert_size_metrics.txt \
      H=${bam_path}PicardQC/insert_size_histogram.pdf \
      M=0.5
      
#Running samtools statistics on BAM indeces
samtools index ${bam_path}${bam_file_name}
samtools idxstats ${bam_path}${bam_file_name} > ${bam_path}PicardQC/indexStats.txt

#Running RNA-seQC

md_bam_file=${bam_path}markedDup/Aligned.sortedByCoord.out.bam
mkdir -p -- ${bam_path}RNA_SeQC



samtools index ${md_bam_file}

javaPath="/home/vahid/.conda/pkgs/java-1.7.0-openjdk-cos6-x86_64-1.7.0.131-h06d78d4_0/x86_64-conda_cos6-linux-gnu/sysroot/usr/lib/jvm/java-1.7.0-openjdk-1.7.0.131.x86_64/jre/bin/"
${javaPath}java -jar /data/vahid/RNAseqPipeLine/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
    -s ${sample_id},${md_bam_file},${sample_id} \
    -t ${refGTF} \
    -r ${refGenome} \
    -noDoC \
    -strictMode \
    -o ${bam_path}RNA_SeQC \
    -gatkFlags --allow_potentially_misencoded_quality_scores \
    -singleEnd no
    

#Running RSEM

mkdir -p -- ${bam_path}RSEM
/data/vahid/RNAseqPipeLine/RSEMpkg/RSEM-1.2.25/rsem-calculate-expression --num-threads 4 \
        --fragment-length-max 1000 \
        --estimate-rspd \
        --no-bam-output \
        --paired-end \
        --bam ${bam_path}Aligned.toTranscriptome.out.bam \
        ${RSEMpath} ${bam_path}RSEM/RSEM

#Running HT-Seq
mkdir -p -- ${bam_path}HTseq
htseq-count -f bam -r pos -s no -t gene ${bam_path}${bam_file_name} ${refGTF} > ${bam_path}HTseq/HTseqGeneUnion.txt &
sleep 1m
htseq-count -f bam -r pos -s no -t transcript ${bam_path}${bam_file_name} ${refGTF} > ${bam_path}HTseq/HTseqTranscriptUnion.txt &
sleep 1m
htseq-count -f bam -r pos -s no -t exon ${bam_path}${bam_file_name} ${refGTF} > ${bam_path}HTseq/HTseqExonUnion.txt &
sleep 1m
htseq-count -f bam -r pos -s no -t gene -m intersection-strict ${bam_path}${bam_file_name} ${refGTF} > ${bam_path}HTseq/HTseqGeneIntersect.txt &
sleep 1m
htseq-count -f bam -r pos -s no -t transcript -m intersection-strict ${bam_path}${bam_file_name} ${refGTF} > ${bam_path}HTseq/HTseqTranscriptIntersect.txt &
sleep 1m
htseq-count -f bam -r pos -s no -t exon -m intersection-strict ${bam_path}${bam_file_name} ${refGTF} > ${bam_path}HTseq/HTseqExonIntersect.txt
sleep 30m



#rm -- ${bam_path}fastq1.gz
#rm -- ${bam_path}fastq2.gz
rm -rf -- ${bam_path}_STARpass1
rm -rf -- ${bam_path}_STARgenome

rm -rf -- ${bam_path}markedDup

rm -- ${bam_path}Aligned.toTranscriptome.out.bam
rm -- ${bam_path}Aligned.sortedByCoord.out.bam
rm -- ${bam_path}Aligned.sortedByCoord.out.bam.bai