from snakemake.remote import FTP


FTP = FTP.RemoteProvider()


configfile: "config.yaml"


rule all:
    input:
        "ref/mus_musculus_38_91/annotation.gtf",
        "ref/mus_musculus_38_91/genome.fa",
        expand("{sample}_{group}.fastq.gz", 
               group=['R1', 'R2'], sample=["sample1", "sample2"])


rule annotation:
    input:
        FTP.remote("ftp.ensembl.org/pub/release-91/gtf/mus_musculus/Mus_musculus.GRCm38.91.gtf.gz", static=True, keep_local=True)
    output:
        "ref/mus_musculus_38_91/annotation.gtf"
    shell:
        "zgrep -P ^{wildcards.chrom} {input} > {output}"


rule genome:
    input:
        FTP.remote("ftp.ensembl.org/pub/release-91/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz", static=True, keep_local=True)
    output:
        "ref/mus_musculus_38_91/genome.fa"
    shell:
        "gzip -d -c {input} > {output}"


rule reads:
    output:
        "{sample}_{group}.fastq.gz"
    params:
        url=config["bam"],
        seed=lambda wildcards: abs(hash(wildcards.sample)) % 10000
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools bam2fq -1 {output[0]} -2 {output[1]} "
        "<(samtools view -b -s{params.seed}.2 {params.url} chr{wildcards.chrom})"
