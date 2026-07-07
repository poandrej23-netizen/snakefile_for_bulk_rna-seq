# snakefile_for_bulk_rna-seq
Snakefile for bulk rna-seq analysis (fastqc, STAR, DeSeq2, Gene Ontology)

Requires snakemake to be used. Supports comparison of different groups of samples. 
Generates table of differentially expressed genes and GO enrichment plots with top GO terms for up- and downregulated genes. 
Tresholds for log2FC and padj, genome version and sample names can be regulated directly in Snakefile
