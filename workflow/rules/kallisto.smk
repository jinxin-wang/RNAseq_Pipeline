## quantifying abundances of transcripts from RNA-Seq data
rule quantification_with_kallisto:
    input:
        reads = ["fastq/{sample}_1.fastq.gz", "fastq/{sample}_2.fastq.gz"] if config["paired"] == True  else "fastq/{sample}_0.fastq.gz",
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
        "  -o {output.dir} {input.reads} 2> {log} "
    
rule concatenate_abundances:
    input:
        expand("kallisto_targetID_count/{sample}/abundance.tsv", sample=SAMPLES),
    output:
        tech   = temp("kallisto_targetID_count/tech.txt"),
        tpm    = temp("kallisto_targetID_count/tpm.txt"),
        header = temp("kallisto_targetID_count/header.tsv"),
        tsv    = temp("kallisto_targetID_count/transcript_tpms_lung_multirigional.tsv")
    log:
        "logs/kallisto_targetID_count/transcript_tpms_lung_multirigional.log"
    threads: 1
    params:
        queue = "shortq"
    resources:
        mem_mb = 20480
    run:
        shell("paste {input} | cut -f 1 > {output.tech}"),
        shell("paste {input} | awk -F'\t' '{for(i=5;i<=NF;i+=5){printf \"%s\t\",$i;} print \"\"}' | paste {output.tech} - > {output.tpm}"),
        shell("ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print \"\t$1\"}' | perl -ne 'print \"target_id$_\n\"' > {output.header}"),
        shell("cat {output.header} {output.tpm} | grep -w -v \"tpm\" | sed 's/_R1.fastq//g' | sed 's/[[:space:]]*$//' > {output.tsv}")
