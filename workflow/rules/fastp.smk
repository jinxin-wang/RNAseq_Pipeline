## A rule to clean fastq files with fastp
rule fastp:
    input:
        fastq_0="fastq/{sample}_0.fastq.gz" if config["paired"] == False, 
        fastq_1="fastq/{sample}_1.fastq.gz" if config["paired"] == True , 
        fastq_2="fastq/{sample}_2.fastq.gz" if config["paired"] == True ,
    output:
        fastq_clean_0='fastq_clean/{sample}_0.fastq.gz' if config["paired"] == False,
        fastq_clean_1='fastq_clean/{sample}_1.fastq.gz' if config["paired"] == True ,
        fastq_clean_2='fastq_clean/{sample}_2.fastq.gz' if config["paired"] == True ,
        html_report='fastp_reports/{sample}_fastp_report.html',
        json_report='fastp_reports/{sample}_fastp_report.json'
    log:
        "logs/fastp/{sample}_fastp.html.log"
    params:
        queue = "mediumq",
        fastp = config["fastp"]["app"],
        adapters = config["fastp"]["adapters"]
    threads : 16
    resources:
        mem_mb = 51200
    run:
        if config["paired"] == False :
            print("[Message: rule fastp] Clean fastq files with fastp single read.")
            shell('{params.fastp} --thread {threads} --dont_overwrite -i {input.fastq_0} -o {output.fastq_clean} --compression 9 --adapter_fasta {params.adapters} --trim_poly_g --trim_poly_x --length_required 25 --overrepresentation_analysis  --html {output.html_report} --json {output.json_report} 2> {log}')
        elif config["paired"] == True :
            print("[Message: rule fastp] Clean fastq files with fastp paired-end.")
            shell('{params.fastp} --thread {threads} --dont_overwrite -i {input.fastq_1} -o {output.fastq_clean_1} -I {input.fastq_2} -O {output.fastq_clean_2} --compression 9 --adapter_fasta {params.adapters} --trim_poly_g --trim_poly_x --length_required 25 --overrepresentation_analysis  --html {output.html_report} --json {output.json_report} 2> {log}')
        else:
            print("[Message: rule fastp] Unable to recognize the configuration value of 'paired': %s"%(str(config["paired"])))
