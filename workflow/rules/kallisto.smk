## quantifying abundances of transcripts from RNA-Seq data
rule quantification_with_kallisto:
    input:
        reads = "fastq/{sample}_1.fastq.gz", "fastq/{sample}_2.fastq.gz" if config["paired"] else "fastq/{sample}_0.fastq.gz"
    output:
        dir   = "kallisto_targetID_count/{sample}/",
        h5    = "kallisto_targetID_count/{sample}/abundance.h5",
        table = "kallisto_targetID_count/{sample}/abundance.tsv",
    log:
        "logs/kallisto_targetID_count/{sample}.log",
    threads: 16
    params:
        queue = "mediumq",
        kallisto = config["kallisto"]["app"],
        index = config["kallisto"]["humain"]["index"] if config["samples"] == "humain" else config["kallisto"]["mouse"]["index"]
    resources:
        mem_mb = 51200
    shell:
        "{params.kallisto} quant --verbose -t {threads} "
        "  -i {params.index} "
        "  -o {output.dir} {input.reads} "
    
