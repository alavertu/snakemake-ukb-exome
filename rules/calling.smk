
rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"]
    output:
        gvcf=protected("data/called/{sample}.g.vcf.gz")
    log:
        "data/logs/gatk/haplotypecaller/{sample}.log"
    params:
        java_opts="-Xmx50g",
        extra=" -L " + config["processing"].get("restrict-regions") + " --interval-padding " + config["processing"].get("region-padding")
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
        gvcfs=expand("data/called/{sample}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="data/called/all.g.vcf.gz"
    log:
        "data/logs/gatk/combinegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="data/called/all.g.vcf.gz"
    output:
        vcf="data/genotyped/all.vcf.gz"
    params:
          java_opts="-Xmx50g",
          extra=" -L " + config["processing"].get("restrict-regions") + " --interval-padding " + config["processing"].get("region-padding")
    log:
        "data/logs/gatk/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"
