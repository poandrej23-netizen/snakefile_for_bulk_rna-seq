import os

GENOME_DIR = "/mnt/storage/lab3/Abramov/reference/hg38"
GTF_FILE   = "/mnt/storage/lab3/Abramov/reference/hg38/gencode.v46.basic.annotation.gtf"
THREADS    = 16

GLOBAL_TMPDIR = os.path.abspath("./snakemake_tmp")
os.makedirs(GLOBAL_TMPDIR, exist_ok=True)
os.environ["TMPDIR"] = GLOBAL_TMPDIR

COMPARISONS = {
    "HEK_Trim": {
        "WT": ["HEK_WT1", "HEK_WT2", "HEK_WT3"],
        "KO": ["HEK_Trim_KO1", "HEK_Trim_KO2", "HEK_Trim_KO3"]
    }
}

ALL_SAMPLES = []
for comp in COMPARISONS.values():
    ALL_SAMPLES.extend(comp["WT"])
    ALL_SAMPLES.extend(comp["KO"])
ALL_SAMPLES = list(set(ALL_SAMPLES))
SAMPLES = ALL_SAMPLES

PADJ_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0

rule all:
    input:
        "counts/final_counts.tsv",
        expand("qc/{sample}_fastqc.html", sample=ALL_SAMPLES),
        expand("deseq/{comparison}_processed.done", comparison=COMPARISONS.keys())

rule fastqc:
    input:
        "data/{sample}.fastq.gz"
    output:
        html="qc/{sample}_fastqc.html",
        zip="qc/{sample}_fastqc.zip"
    threads: 4
    shell:
        """
        mkdir -p qc
        export _JAVA_OPTIONS="-Djava.io.tmpdir={GLOBAL_TMPDIR}"
        
        fastqc -t {threads} -o qc {input}
        """

rule star_align:
    input:
        "data/{sample}.fastq.gz"
    output:
        "bams/{sample}Aligned.sortedByCoord.out.bam"
    threads: THREADS
    shell:
        """
        mkdir -p bams
        STAR --runThreadN {threads} \
             --genomeDir {GENOME_DIR} \
             --readFilesCommand zcat \
             --readFilesIn {input} \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix bams/{wildcards.sample} \
             --outFilterMismatchNoverLmax 0.05
        """

rule htseq_count:
    input:
        bam="bams/{sample}Aligned.sortedByCoord.out.bam"
    output:
        "counts/{sample}.counts"
    shell:
        """
        mkdir -p counts
        htseq-count -f bam -r pos -s no {input.bam} {GTF_FILE} > {output}
        """


rule combine_counts:
    input:
        expand("counts/{sample}.counts", sample=SAMPLES)
    output:
        "counts/final_counts.tsv"
    shell:
        """
        TMP_DIR=$(mktemp -d -p .)
       
        cut -f1 {input[0]} > $TMP_DIR/gene_ids.txt

        header="GeneID"
        count_files=""
        for f in {input}; do
            sample=$(basename $f .counts)
            header="$header\\t$sample"
            cut -f2 $f > $TMP_DIR/${{sample}}_counts.txt
            count_files="$count_files $TMP_DIR/${{sample}}_counts.txt"
        done
       
        echo -e "$header" > {output}
        paste $TMP_DIR/gene_ids.txt $count_files >> {output}
       
        rm -rf $TMP_DIR
        """

rule run_deseq:
    input:
        counts = "counts/final_counts.tsv"
    output:
        res_annotated = "deseq/{comparison}_DESeq2_results_annotated.csv",
        go_up = "deseq/{comparison}_GO_enrichment_up.csv",
        go_down = "deseq/{comparison}_GO_enrichment_down.csv",
        plot_up = "deseq/{comparison}_GO_dotplot_up.pdf",
        plot_down = "deseq/{comparison}_GO_dotplot_down.pdf",
        done = touch("deseq/{comparison}_processed.done")
    params:
        wt_samples = lambda wildcards: COMPARISONS[wildcards.comparison]["WT"],
        ko_samples = lambda wildcards: COMPARISONS[wildcards.comparison]["KO"],
        padj = PADJ_CUTOFF,
        log2fc = LOG2FC_CUTOFF,
        cache_dir = GLOBAL_TMPDIR
    conda:
        "envs/deseq2.yaml"
    script:
        "script/deseq.R"
