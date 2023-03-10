#!/bin/bash

/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake -s snakemake.smk -c 'sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M -p {params.queue}' --jobs 20 --rerun-incomplete --config samples="humain" paired=True
