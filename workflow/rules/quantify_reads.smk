## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_transcriptID:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        transcript_id_count = "HTSeq_transcriptID_count/{sample}_transcriptID_count.table"
    log:
        "logs/HTSeq_transcriptID_count/{sample}_transcriptID_count.log"
    threads: 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i transcript_id"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.transcript_id_count} 2> {log}"
    
## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_genesNAME:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        gene_name_count = "HTSeq_geneNAME_count/{sample}_geneNAME_count.table"
    log:
        "logs/HTSeq_geneNAME_count/{sample}_geneNAME_count.log"
    threads : 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i gene_name"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.gene_name_count} 2> {log}"

## A rule to quantify reads per annotation, with HTSeq
rule quantification_with_HTSeq_transcriptNAME:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        htseq_gtf = config["HTSEQ_GTF"]
    output:
        transcript_name_count = "HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.table"
    log:
        "logs/HTSeq_transcriptNAME_count/{sample}_transcriptNAME_count.log"
    threads : 1
    params:
        queue = "mediumq",
        htseq_count = config["APP_HTSEQ_COUNT"],
        strandness  = config["HTSEQ_STRANDNESS"]
    resources:
        mem_mb = 20480
    shell:
        "{params.htseq_count}  "
        " -f bam"
        " -i transcript_name"
        " -r pos"
        " {params.strandness}"
        " {input.bam}"
        " {input.htseq_gtf}"
        " > {output.transcript_name_count} 2> {log}"
