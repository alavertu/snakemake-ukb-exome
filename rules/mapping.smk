rule revert_cram_to_sam:
    input:
         unpack(get_cram)
    output:
          "sam_files/{sample}-{unit}.sam"
    params:
          "-h -t /oak/stanford/groups/rbaltman/references/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" # optional params string
    wrapper:
           "0.43.1/bio/samtools/view"

rule sort_sam:
    input:
         "sam_files/{sample}-{unit}.sam"
    output:
          "sam_files/{sample}-{unit}.sorted.sam"
    shell:
        "samtools collate -o {output} {input}"

rule sam_to_fastq:
    input:
         "sam_files/{sample}-{unit}.sorted.sam"
    output:
          single="fastq/{sample}-{unit}.sorted.single",
          pe1="fastq/{sample}-{unit}.sorted.1.fq",
          pe2="fastq/{sample}-{unit}.sorted.2.fq"
    shell:
         "samtools fastq -s {output.single} -1 {output.pe1} -2 {output.pe2} {input}"

rule trim_reads_pe:
    input:
        unpack(get_fastq)
    output:
        r1=temp("trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("trimmed/{sample}-{unit}.2.fastq.gz"),
        trimlog="trimmed/{sample}-{unit}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        **config["params"]["trimmomatic"]["pe"]
    log:
        "logs/trimmomatic/{sample}-{unit}.log"
    wrapper:
        "0.30.0/bio/trimmomatic/pe"


rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        temp("mapped/{sample}-{unit}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}-{unit}.log"
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
        "mapped/{sample}-{unit}.sorted.bam"
    output:
        bam=temp("dedup/{sample}-{unit}.bam"),
        metrics="qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
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
        bam=protected("recal/{sample}-{unit}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"
