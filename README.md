## [Institut Gustave Roussy] RNAseq analysis pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

## Usage on flamingo

- Step 0. clone github workflow project on flamingo
```
$ ssh username@flamingo.intra.igr.fr
$ git clone https://github.com/jinxin-wang/RNAseq_Pipeline.git
```
- Step 1. create conda envirements 
```
$ conda env create -f RNA_env.txt
$ conda env create -f pipeline_GATK_2.1.4_conda_env.txt
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
$ ./run.sh
