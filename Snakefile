include: "rules/common.smk"

##### Target rules #####

rule all:
    input:
        "data/annotated/all.vcf.gz",
        "data/qc/multiqc.html",
        "data/plots/depths.svg",
        "data/plots/allele-freqs.svg"


##### Modules #####

include: "rules/mapping.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
