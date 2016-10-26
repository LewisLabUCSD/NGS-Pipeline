# NewPipeline
Pipelines to process NGS or Pacbio data
---------------------------------------

## Method to run these pipelines.
1. define all parameters in the corresponding parameter file in parameters folder.
2. In bash terminal, run the followsing command
	nohup python pipeline.py parameter.yaml > log.txt &
3. Press enter

* Finished Pipeline
	* RNAseq_count: quantify number of reads mapping to each gene
	* GATK_RNA_CHO: call variants for RNAseq
	* SV_Pacbio_PBHoney: call structure variation for Pacbio data using PBHoney
	* SV_Pacbio_Sniffle: call structure variation for Pacbio data using Sniffle
	

	