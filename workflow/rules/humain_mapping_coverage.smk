## A rule to check mapping coverage with rseqc
rule rseqc_geneBody_coverage:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
        rseqc_ref = config["RSEQC_REF"]
    output:
        o1 = "rseqc_geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf"
        # o2 = "rseqc_geneBody_coverage/{sample}.geneBodyCoverage.r",
        # o3 = "rseqc_geneBody_coverage/{sample}.geneBodyCoverage.txt"
    log:
        "logs/rseqc_geneBody_coverage/{sample}.log"
    threads : 1
    params:
        queue = "mediumq",
        genbd = config["APP_GENE_BD_CVR"]
    resources:
        mem_mb = 5120
    shell:
        "cd rseqc_geneBody_coverage/ ;  {params.genbd} -r {input.rseqc_ref} -l 500 -i ../{input.bam} -o {wildcards.sample} 2> ../{log}"
