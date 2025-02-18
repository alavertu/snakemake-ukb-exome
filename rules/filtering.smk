def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="data/genotyped/" + config["run"]["group"] + ".vcf.gz"
    output:
        vcf=temp("data/filtered/" + config["run"]["group"] + ".{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "data/logs/gatk/selectvariants/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/selectvariants"


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="data/filtered/" + config["run"]["group"] + ".{vartype}.vcf.gz"
    output:
        vcf=temp("data/filtered/" + config["run"]["group"] + ".{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        "data/logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf="data/filtered/" + config["run"]["group"] + ".{vartype}.vcf.gz"
    output:
        vcf=temp("data/filtered/" + config["run"]["group"] + ".{vartype}.recalibrated.vcf.gz")
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        "data/logs/gatk/variantrecalibrator/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantrecalibrator"


rule merge_calls:
    input:
        vcf=expand("data/filtered/" + config["run"]["group"] + ".{vartype}.{filtertype}.vcf.gz",
                   vartype=["snvs", "indels"],
                   filtertype="recalibrated"
                              if config["filtering"]["vqsr"]
                              else "hardfiltered")
    output:
        vcf="data/filtered/" + config["run"]["group"] + ".vcf.gz"
    log:
        "data/logs/picard/merge-filtered.log"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
