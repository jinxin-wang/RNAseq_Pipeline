## A rule to check mapping metrics with samtools flagstat
rule samtools_flagstat:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai"
    output:
        "mapping_QC/{sample}_flagstat.txt"
    log:
        "logs/mapping_QC/{sample}.flagstat.log"
    threads : 1
    params:
        queue = "shortq",
        samtools = config["samtools"]["app"],
    resources:
        mem_mb = 2000
    shell:
        "{params.samtools} flagstat {input.bam} > {output} 2> {log}"

