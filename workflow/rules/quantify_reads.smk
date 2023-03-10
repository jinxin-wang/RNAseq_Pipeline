## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_transcriptID:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
    output:
        transcript_id_count = "HTSeq_transcriptID_count/{sample}_transcriptID_count.table"
    log:
        "logs/HTSeq_transcriptID_count/{sample}_transcriptID_count.log"
    threads: 1
    params:
        queue = "mediumq",
        htseq = config["htseq"]["app"],
        strandness = config["htseq"]["strandness"],
        gtf = config["htseq"]["humain"]["gtf"] if config["samples"] == "humain" else config["htseq"]["mouse"]["gtf"] ,
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq} "
        " -f bam"
        " -i transcript_id"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {params.gtf}"
        " > {output.transcript_id_count} 2> {log}"
    
## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_genesNAME:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
    output:
        gene_name_count = "HTSeq_geneNAME_count/{sample}_geneNAME_count.table"
    log:
        "logs/HTSeq_geneNAME_count/{sample}_geneNAME_count.log"
    threads : 1
    params:
        queue = "mediumq",
        htseq = config["htseq"]["app"],
        strandness = config["htseq"]["strandness"],
        gtf = config["htseq"]["humain"]["gtf"] if config["samples"] == "humain" else config["htseq"]["mouse"]["gtf"] ,
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq}  "
        " -f bam"
        " -i gene_name"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {params.gtf}"
        " > {output.gene_name_count} 2> {log}"

## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_transcriptNAME:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
    output:
        transcript_name_count = "HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.table"
    log:
        "logs/HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.log"
    threads : 1
    params:
        queue = "mediumq",
        htseq = config["htseq"]["app"],
        strandness = config["htseq"]["strandness"],
        gtf = config["htseq"]["humain"]["gtf"] if config["samples"] == "humain" else config["htseq"]["mouse"]["gtf"] ,
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq}  "
        " -f bam"
        " -i transcript_name"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {params.gtf}"
        " > {output.transcript_name_count} 2> {log}"
