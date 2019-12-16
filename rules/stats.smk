rule vcf_to_tsv:
    input:
        "data/annotated/all.vcf.gz"
    output:
        report("data/tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


rule plot_stats:
    input:
        "data/tables/calls.tsv.gz"
    output:
        depths=report("data/plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("data/plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/plot-depths.py"
