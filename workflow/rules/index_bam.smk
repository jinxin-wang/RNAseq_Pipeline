## A rule to generate bam index with samtools
rule indexbam:
    input:
        bam = "bam/{sample}.bam"
    output:
        bai = "bam/{sample}.bam.bai"
    threads : 1
    resources:
        mem_mb = 2000
    params:
        queue = "shortq"
    log:
        "logs/bam/{sample}.bam.bai.log"
    shell:
        "samtools index {input} 2> {log}"
