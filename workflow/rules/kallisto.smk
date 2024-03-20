## quantifying abundances of transcripts from RNA-Seq data
rule quantification_with_kallisto:
    input:
        reads = ["fastq/{sample}_1.fastq.gz", "fastq/{sample}_2.fastq.gz"] if config["paired"] == True  else "fastq/{sample}_0.fastq.gz",
    output:
        dir   = directory("kallisto_targetID_count/{sample}/"),
        h5    = "kallisto_targetID_count/{sample}/abundance.h5",
        table = "kallisto_targetID_count/{sample}/abundance.tsv",
    log:
        "logs/kallisto_targetID_count/{sample}.log",
    threads: 4
    params:
        queue = "shortq",
        kallisto = config["kallisto"]["app"],
        index = config["kallisto"][config["samples"]]["index"], # if config["samples"] == "humain" else config["kallisto"]["mouse"]["index"]
    resources:
        mem_mb = 10240
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
        trans  = "kallisto_targetID_count/tpm_transcript.tsv",
    log:
        "logs/kallisto_targetID_count/concatenate_abundances.log",
    threads: 1,
    params:
        queue = "shortq",
    resources:
        mem_mb = 2048
    shell:
        "paste {input} | cut -f 1 > {output.tech} ; "
        "paste {input} | awk -F'\t' '{{for(i=5;i<=NF;i+=5){{printf \"%s\t\",$i;}} print \"\"}}' | paste {output.tech} - > {output.tpm} ; "
       	"grep -l \"ENST\" {input} | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){{print \"\t$1\"}}' | perl -ne 'print \"target_id$_\n\"' > {output.header} ; "
        "cat {output.header} {output.tpm} | grep -w -v \"tpm\" | sed 's/_R1.fastq//g' | sed 's/[[:space:]]*$//' > {output.trans} ; "


rule build_genes_tables:
    input:
        tech   = "kallisto_targetID_count/tech.txt",
        tpm    = "kallisto_targetID_count/tpm.txt",
        header = "kallisto_targetID_count/header.tsv",
        trans  = "kallisto_targetID_count/tpm_transcript.tsv",
    output:
        genes  = "kallisto_targetID_count/tpm_genes.tsv",
        kassandra = "kallisto_targetID_count/tpm_kassandra.tsv",
    log:
        "logs/kallisto_targetID_count/build_genes_tables.log",
    threads: 1,
    params:
        queue = "shortq",
        trans2genes = config["kallisto"]["trans2genes"],
        transfilter = config["kallisto"]["transfilter"],
        id2genes    = config["kallisto"][config["samples"]]["id2genes"],
        genes_expr  = config["kallisto"][config["samples"]]["genes_expr"],
    resources:
        mem_mb = 2048
    conda: 'routine'
    shell:
        "echo 'starting {params.trans2genes}' >> {log} ;"
        "Rscript {params.trans2genes} {input.trans} {params.id2genes} {output.genes} 2>> {log} ; "
        "echo 'starting {params.transfilter}' >> {log} ;"
        "Rscript {params.transfilter} {input.trans} {params.genes_expr} {params.id2genes} {output.kassandra} 2>> {log} ; "

