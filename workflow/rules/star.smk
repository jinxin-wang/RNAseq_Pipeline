## A rule to map single-end RNA sample using STAR
if config["Paired"] == False: 
    rule star_single:
        input:
            reads1 = "fastq_clean/{sample}_0.fastq.gz"
        output:
            bam = "bam/{sample}.bam"
        params:
            queue = "mediumq",
            star  = config["APP_STAR"],
            index = config["STAR_INDEX"],
            tmp_directory_for_index = "TMP_STAR/tmp_genome_for_star_{sample}/",
            genome_fasta = config["GENOME_FASTA"],
            star_sjdbOverhang = config["STAR_SJDBOVERHANG"],
            star_gtf = config["STAR_GTF"]
        log:
            "logs/bam/{sample}.bam.log"
        threads: 16
        resources:
            mem_mb = 51200
        shell:
            "mkdir -p {params.tmp_directory_for_index};"
            "{params.star}"
            " --readFilesCommand zcat"
            " --outStd Log"
            " --runThreadN {threads}"
            " --genomeDir {params.index}"
            " --readFilesIn {input.reads1}"
            " --sjdbOverhang {params.star_sjdbOverhang}"
            " --sjdbGTFfile {params.star_gtf}"
            " --twopassMode Basic"
            " --outSAMtype BAM SortedByCoordinate"
            " --outSAMunmapped Within"
            " --outSAMstrandField intronMotif"
            " --outSAMmultNmax 5"
            " --outSJfilterReads Unique"
            " --quantMode GeneCounts"
            " --outFileNamePrefix {params.tmp_directory_for_index}{wildcards.sample}_TMP_; "
            " mv {params.tmp_directory_for_index}{wildcards.sample}_TMP_Aligned.sortedByCoord.out.bam {output.bam} 2>> {log};"

## A rule to map paired-end RNA sample using STAR
## [J. WANG] ATTENTION :
##      if threads number is more than 12, then you need to add the command "ulimit -n 10000" 
##      at the begin of the shell script. 
if config["Paired"] == True:
    rule star_paired:
        input:
            reads1 = "fastq_clean/{sample}_1.fastq.gz",
            reads2 = "fastq_clean/{sample}_2.fastq.gz"
        output:
            bam = "bam/{sample}.bam"
        params:
            queue = "mediumq",
            star  = config["APP_STAR"],
            index = config["STAR_INDEX"],
            tmp_directory_for_index = "TMP_STAR/tmp_genome_for_star_{sample}/",
            genome_fasta = config["GENOME_FASTA"],
            star_sjdbOverhang = config["STAR_SJDBOVERHANG"],
            star_gtf = config["STAR_GTF"],
            start_outSAMmultNmax = config["START_OUTSAMMULTNMAX"]
        log:
            "logs/bam/{sample}.bam.log"
        threads: 32
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
            " --readFilesIn {input.reads1} {input.reads2}"
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
 
