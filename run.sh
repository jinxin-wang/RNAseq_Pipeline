#!/bin/bash

set -e ; 
source ~/.bashrc ;

conda activate routine ; 

/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake -c 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M -p {params.queue}' --jobs 20 --rerun-incomplete --config samples="human" paired=True do_qc=False do_fastp=False do_bam=False do_mqc=False do_htseq=False do_rseqc=False do_kallisto=True --nt
