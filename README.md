## [Institut Gustave Roussy] RNAseq analysis pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-=5.23.0-brightgreen.svg)](https://snakemake.github.io)

## Usage on flamingo

- Step 0. clone github workflow project on flamingo
```
$ ssh username@flamingo.intra.igr.fr
$ cd /workflowpath
$ git clone https://github.com/jinxin-wang/RNAseq_Pipeline.git
```
- Step 0. to update the code
```
$ cd /workflowpath/RNAseq_Pipeline
$ git pull
```
- Step 1. create conda envirements 
```
$ cd RNAseq_Pipeline
$ conda env create --file envs/RNA_env.txt --name RNA
$ conda env create --file envs/pipeline_GATK_2.1.4_conda_env.txt --name pipeline_GATK_2.1.4_V2
```
- Step 2. deploy workflow
```
$ cd /mnt/beegfs/scratch/username/yourprojectdir
$ mkdir projectname
$ cd projectname
$ ln -s /workflowpath/RNAseq_Pipeline/workflow .
$ cp /workflowpath/RNAseq_Pipeline/run.sh .
```
- Step 3. configure workflow
```
$ emacs -nw run.sh
```
- Step 4. run workflow
```
$ conda activate pipeline_GATK_2.1.4_V2
$ ./run.sh
