# NewPipeline
Pipelines to process NGS or Pacbio data
---------------------------------------
## Softwares needed
All the softwares needed for each pipeline is listed in the corresponding parameter files.
## Method to run these pipelines.
#### All fastq files should be compressed. Paired end files should end with _1.fq.gz, _2.fq.gz or _1.fastq.gz,_2.fastq.fz. Single end files should end with _1.fq.gz
#### gff file and refernce fa file should not be zipped.
#### STAR takes a lot of memory(30-50 GB) each run, so don't run more than 2 STAR in parallel at each batch.
1. define all parameters in the corresponding parameter file in parameters folder.
2. In bash terminal, run the followsing command:
	* nohup python /path/to/pipeline.py /path/to/parameter.yaml > log.txt & (don't forget to include & symbol)
3. Press enter

* Finished Pipeline
	* RNAseq_count: quantify number of reads mapping to each gene
	* GATK_RNA_CHO: call variants for RNAseq
	* SV_Pacbio_PBHoney: call structure variation for Pacbio data using PBHoney
	* SV_Pacbio_Sniffle: call structure variation for Pacbio data using Sniffle
	

	
