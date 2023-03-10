## A rule to map single-end or paired-end RNA sample using STAR
## [J. WANG] ATTENTION :
##      if threads number is more than 12, then you need to add the command "ulimit -n 10000" 
##      at the begin of the shell script. 
rule star:
    input:
        reads = "fastq_clean/{sample}_1.fastq.gz", "fastq_clean/{sample}_2.fastq.gz" config["paired"] == True else "fastq_clean/{sample}_0.fastq.gz"
    output:
        bam = "bam/{sample}.bam"
    params:
        queue = "mediumq",
        tmp_directory_for_index = "TMP_STAR/tmp_genome_for_star_{sample}/",
        star  = config["star"]["app"],
        index = config["star"]["humain"]["index"] if config["samples"] = "humain" else config["star"]["mouse"]["index"],
        genome_fasta = config["star"]["humain"]["genome_fasta"] if config["samples"] = "humain" else config["star"]["mouse"]["genome_fasta"],
        star_sjdbOverhang = config["star"]["humain"]["sjdboverhang"] if config["samples"] = "humain" else config["star"]["mouse"]["sjdboverhang"],
        star_gtf = config["star"]["humain"]["gtf"] if config["samples"] = "humain" else config["star"]["mouse"]["gtf"],
    log:
        "logs/bam/{sample}.bam.log"
    threads: 24
    resources:
        mem_mb = 102400
    shell:
        "ulimit -n 10000 2>> {log}; "
        "mkdir -p {params.tmp_directory_for_index} 2>> {log} ; "
        "{params.star}"
        " --readFilesCommand zcat"
        " --outStd Log"
        " --runThreadN {threads}"
        " --genomeDir {params.index}"
        " --readFilesIn {input.reads}"
        " --sjdbOverhang {params.star_sjdbOverhang}"
        " --sjdbGTFfile {params.star_gtf}"
        " --twopassMode Basic"
        " --outSAMtype BAM SortedByCoordinate"
        " --outSAMunmapped Within"
        " --outSJfilterReads Unique"
        " --quantMode GeneCounts"
        " --outSAMstrandField intronMotif"
        " --outSAMmultNmax 5"
        " --outFileNamePrefix {params.tmp_directory_for_index}{wildcards.sample}_TMP_ 2>> {log} ; "
        " mv {params.tmp_directory_for_index}{wildcards.sample}_TMP_Aligned.sortedByCoord.out.bam {output.bam} 2>> {log};"
 
