FROM continuumio/miniconda
MAINTAINER Shangzhong Li (shl198@eng.ucsd.edu)

CMD ["date"]
ENTRYPOINT ["python", "/root/Code/NGS-Pipeline/Salmon_quant.py"]


RUN mkdir ~/Code && cd ~/Code && \
	git clone --single-branch -b ngs_b1 https://github.com/LewisLabUCSD/NGS-Pipeline.git && \
	apt-get update && conda install -c bioconda ruffus && \
	conda install natsort && \
	conda install -c conda-forge sarge && \
	conda install -c bioconda trimmomatic && \
	conda install -c bioconda fastqc && \
	conda install -c bioconda salmon && \
	conda install -c conda-forge pandas && \
	conda install -c conda-forge pyyaml 
	







