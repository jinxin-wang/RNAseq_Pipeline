## A rule to generate fastq quality control
rule fastqc:
    input:
        fastq='fastq/{fastq_sample}.fastq.gz',
    output:
        'fastq_QC/{fastq_sample}_fastqc.html',
        'fastq_QC/{fastq_sample}_fastqc.zip'
    log:
        "logs/fastq_QC/{fastq_sample}_fastqc.html.log"
    threads : 4
    resources:
        mem_mb = 51200
    params:
        queue = "shortq",
        fastqc = config["fastqc"]["app"],
        adapters = config["fastqc"]["adapters"],
    shell:
        'module load java ; '
        '{params.fastqc} -t {threads} -a {params.adapters} -o fastq_QC/ {input.fastq} 2> {log}'


## A rule to generate fastq quality control for cleaned fastq
rule fastqc_clean:
    input:
        fastq='fastq_clean/{fastq_sample}.fastq.gz'
    output:
        'fastq_QC_clean/{fastq_sample}_fastqc.html',
        'fastq_QC_clean/{fastq_sample}_fastqc.zip'
    log:
        "logs/fastq_QC_clean/{fastq_sample}_fastqc.html.log"
    threads : 4
    resources:
        mem_mb = 51200
    params:
        queue = "shortq",
        fastqc = config["fastqc"]["app"],
        adapters = config["fastqc"]["adapters"]
    shell:
        'module load java ; '
        '{params.fastqc} -t {threads} -a {params.adapters} -o fastq_QC_clean/ {input} 2> {log}'
