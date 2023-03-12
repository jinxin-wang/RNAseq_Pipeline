## A rule to check transcript integrity number with rseqc 
rule rseqc_tin:
     input:
         bam = "bam/{sample}.bam",
         bai = "bam/{sample}.bam.bai",
     output:
         o1 = "rseqc_tin/{sample}.tin.xls",
         o2 = "rseqc_tin/{sample}.summary.txt"
     log:
         "logs/rseqc_tin/{sample}.log"
     threads : 1
     params:
         queue = "mediumq",
         tin = config["tin"]["app"],
         ref = config["tin"]["humain"]["ref"],
     resources:
         mem_mb = 5120
     shell:
         "cd rseqc_tin/ ; {params.tin} -r {params.ref} -c 30 -i ../{input.bam} 2> ../{log}"
        
## A rule to check transcript integrity number with rseqc
rule rseqc_readDuplication:
    input:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai",
    output:
        o1 ="rseqc_read_duplication/{sample}.DupRate_plot.pdf",
        o2 ="rseqc_read_duplication/{sample}.pos.DupRate.xls",
        o3 ="rseqc_read_duplication/{sample}.seq.DupRate.xls",
        o4 ="rseqc_read_duplication/{sample}.DupRate_plot.r"
    log:
        "logs/rseqc_readDuplication/{sample}.log"
    threads : 1
    params:
        queue = "mediumq",
        read_dup = config["read_dup"]["app"]
    resources:
        mem_mb = 51200
    shell:
        "cd rseqc_read_duplication/; {params.read_dup} -i ../{input.bam} -o {wildcards.sample} 2> ../{log}"
