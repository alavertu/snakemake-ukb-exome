include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "data/genotyped/all.vcf.gz",

##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/qc.smk"
