rule revert_cram_to_sam:
    input:
         "/scratch/groups/rbaltman/ukbiobank/aa/{sample}.cram"
    output:
          "data/sam_files/{sample}.sam"
    params:
          "-h -t /oak/stanford/groups/rbaltman/references/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" # optional params string
    wrapper:
           "0.43.1/bio/samtools/view"

rule sort_sam:
    input:
         "data/sam_files/{sample}.sam"
    output:
         "sam_files/{sample}.sorted.sam"
    shell:
        "samtools collate -o {output} {input}"

rule sam_to_fastq:
    input:
         "sam_files/{sample}.sorted.sam"
    output:
          single="fastq/{sample}.sorted.single",
          pe1="fastq/{sample}.sorted.1.fq",
          pe2="fastq/{sample}.sorted.2.fq"
    shell:
         "samtools fastq -s {output.single} -1 {output.pe1} -2 {output.pe2} {input}"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}.1.fastq.gz"),
        r2=temp("trimmed/{sample}.2.fastq.gz"),
        trimlog="trimmed/{sample}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        "logs/trimmomatic/{sample}.log"
    wrapper:
        "0.30.0/bio/trimmomatic/pe"


rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        temp("mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.27.1/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        bam=temp("dedup/{sample}.bam"),
        metrics="qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
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
        bam=protected("recal/{sample}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}.log"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"
