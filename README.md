# NewPipeline
Pipelines to process NGS or Pacbio data
---------------------------------------

## Method to run these pipelines.
#### Paired end files should end with _1.fq.gz, _2.fq.gz or _1.fastq.gz,_2.fastq.fz. Single end files should end with _1.fq.gz
#### STAR takes a lot of memory(30-50 GB) each run, so don't run more than 2 STAR in parallel at each batch.
1. define all parameters in the corresponding parameter file in parameters folder.
2. In bash terminal, run the followsing command:
	* nohup python pipeline.py parameter.yaml > log.txt &
	Or if you are running in screen, try the following command:
	* python pipeline.py parameter.yaml 2>&1 | tee log.txt
3. Press enter

* Finished Pipeline
	* RNAseq_count: quantify number of reads mapping to each gene
	* GATK_RNA_CHO: call variants for RNAseq
	* SV_Pacbio_PBHoney: call structure variation for Pacbio data using PBHoney
	* SV_Pacbio_Sniffle: call structure variation for Pacbio data using Sniffle
	

## Pipeline specific notes

### RNAseq_STARpipeline.sh

* Be aware, this is a bash pipeline and does not manage flow and reruns like ruffus.
	
