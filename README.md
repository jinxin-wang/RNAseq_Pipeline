## [Institut Gustave Roussy] RNAseq analysis pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

## Usage on flamingo

- Step 0. clone github workflow project on flamingo
```
$ ssh username@flamingo.intra.igr.fr
$ cd /workflowpath
$ git clone https://github.com/jinxin-wang/RNAseq_Pipeline.git
$ git pull
```
- Step 1. create conda envirements 
```
$ cd RNAseq_Pipeline
$ conda env create --file envs/RNA_env.txt --name RNAseq
$ conda env create --file envs/pipeline_GATK_2.1.4_conda_env.txt --name pipeline_GATK_2.1.4_V2
```
- Step 2. deploy workflow
```
$ cd /mnt/beegfs/scratch/username/yourprojectdir
$ mkdir projectname
$ cd projectname
$ ln -s /workflowpath/workflow .
$ cp /workflowpath/run.sh .
```
- Step 3. configure workflow
```
$ emacs -nw run.sh
```
- Step 4. run workflow
```
$ conda activate RNAseq
$ ./run.sh
