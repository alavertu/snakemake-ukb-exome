rule revert_cram_to_sam:
    input:
        get_cram
    output:
        temp("data/sam_files/{sample}.sam")
    params:
        "-h -t " + config["ref"]["genome"] + ".fai" # optional params string
    wrapper:
        "0.43.1/bio/samtools/view"

rule sort_sam:
    input:
         "data/sam_files/{sample}.sam"
    output:
         temp("data/sam_files/{sample}.sorted.bam")
    params:
        prefix='data/sam_files/{sample}.sorted'
    shell:
        "samtools collate {input} {params.prefix}"

rule sam_to_fastq:
    input:
         "data/sam_files/{sample}.sorted.bam"
    output:
          single=temp("data/fastq/{sample}.sorted.single"),
          pe1=temp("data/fastq/{sample}.sorted.1.fq"),
          pe2=temp("data/fastq/{sample}.sorted.2.fq")
    shell:
         "samtools fastq -s {output.single} -1 {output.pe1} -2 {output.pe2} {input}"

# rule trim_reads_pe:
#     input:
#         unpack(get_fastq)
#     output:
#         r1=temp("data/trimmed/{sample}.1.fastq.gz"),
#         r2=temp("data/trimmed/{sample}.2.fastq.gz"),
#         r1_unpaired=temp("data/trimmed/{sample}.1.unpaired.fastq.gz"),
#         r2_unpaired=temp("data/trimmed/{sample}.2.unpaired.fastq.gz"),
#         trimlog="data/trimmed/{sample}.trimlog.txt"
#     params:
#         extra=lambda w, output: "-trimlog {}".format(output.trimlog),
#         **config["params"]["trimmomatic"]["pe"]
#     log:
#         "data/logs/trimmomatic/{sample}.log"
#     wrapper:
#         "0.30.0/bio/trimmomatic/pe"


rule map_reads:
    input:
        reads=get_fastq
    output:
        temp("data/mapped/{sample}.bam")
    log:
        "data/logs/bwa_mem/{sample}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        sort="none",
        sort_order="coordinate"
    threads: 4
    wrapper:
        "0.27.1/bio/bwa/mem"

rule sort_aligned_bam:
    input:
         "data/mapped/{sample}.bam"
    output:
          temp("data/mapped/{sample}.sorted.bam")
    params:
        "-m 12G"
    threads:
        4
    wrapper:
           "0.45.1/bio/samtools/sort"


rule mark_duplicates:
    input:
        "data/mapped/{sample}.sorted.bam"
    output:
        bam=temp("data/dedup/{sample}.bam"),
        metrics="data/qc/dedup/{sample}.metrics.txt"
    log:
        "data/logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.26.1/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"]
    output:
        bam=protected("data/recal/{sample}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "data/logs/gatk/bqsr/{sample}.log"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "data/dedup/{sample}.bam"
    output:
        temp("data/dedup/{sample}.bam.bai")
    log:
       "data/logs/picard/dedup/{sample}_indexing.log"
    shell:
        "samtools index {input}"
