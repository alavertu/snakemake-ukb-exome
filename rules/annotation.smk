rule snpeff:
    input:
        "data/filtered/all.vcf.gz",
    output:
        vcf=report("data/annotated/all.vcf.gz", caption="../data/report/vcf.rst", category="Calls"),
        csvstats="data/snpeff/all.csv"
    log:
        "data/logs/snpeff.log"
    params:
        reference=config["ref"]["name"],
        extra="-Xmx6g"
    wrapper:
        "0.27.1/bio/snpeff"
