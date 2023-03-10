## A rule to generate bam index with samtools
rule indexbam:
    input:
        bam = "bam/{sample}.bam"
    output:
        bai = "bam/{sample}.bam.bai"
    threads : 1
    resources:
        mem_mb = 2048
    params:
        queue = "shortq",
        samtools = config["samtools"]["app"],
    log:
        "logs/bam/{sample}.bam.bai.log"
    shell:
        "{params.samtools} index {input} 2> {log}"
