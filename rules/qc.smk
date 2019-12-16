rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="data/qc/fastqc/{sample}.html",
        zip="data/qc/fastqc/{sample}.zip"
    wrapper:
        "0.27.1/bio/fastqc"


rule samtools_stats:
    input:
        "data/recal/{sample}.bam"
    output:
        "data/qc/samtools-stats/{sample}.txt"
    log:
        "data/logs/samtools-stats/{sample}.log"
    wrapper:
        "0.27.1/bio/samtools/stats"


rule multiqc:
    input:
         expand(["data/qc/samtools-stats/{u.sample}.txt",
                 "data/qc/fastqc/{u.sample}.zip",
                 "data/qc/dedup/{u.sample}.metrics.txt"],
                u=units.itertuples()),
         "data/snpeff/all.csv"
    output:
          report("data/qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
       "data/logs/multiqc.log"
    wrapper:
           "0.27.1/bio/multiqc"
