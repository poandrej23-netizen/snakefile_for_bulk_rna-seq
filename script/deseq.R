library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)

cache_dir <- snakemake@params[["cache_dir"]]
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER_CACHE_DIR = cache_dir)
Sys.setenv(HOME = cache_dir)

wt_samples      <- snakemake@params[["wt_samples"]]
ko_samples      <- snakemake@params[["ko_samples"]]
comparison_name <- snakemake@wildcards[["comparison"]]
padj_cutoff     <- snakemake@params[["padj"]]
log2fc_cutoff   <- snakemake@params[["log2fc"]]


counts_data <- read.table(snakemake@input[["counts"]],
                          header = TRUE, row.names = 1,
                          sep = "\t", check.names = FALSE)


samples_to_use <- c(wt_samples, ko_samples)
counts_filtered <- counts_data[, samples_to_use]

condition <- factor(c(rep("WT", length(wt_samples)),
                      rep("KO", length(ko_samples))),
                    levels = c("WT", "KO"))
colData <- data.frame(row.names = samples_to_use, condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData   = colData,
                              design    = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "KO", "WT"))

norm_counts <- counts(dds, normalized = TRUE)
colnames(norm_counts) <- paste0("norm_", colnames(norm_counts))

res_df <- cbind(as.data.frame(res), norm_counts)


res_df$ensembl_id <- gsub("\\.\\d+$", "", rownames(res_df))

write.csv(res_df, snakemake@output[["res_annotated"]], row.names = FALSE)

background_genes <- unique(res_df$ensembl_id)

up_genes <- res_df %>%
    filter(!is.na(padj),
           padj < padj_cutoff,
           log2FoldChange > log2fc_cutoff) %>%
    pull(ensembl_id) %>%
    unique()

down_genes <- res_df %>%
    filter(!is.na(padj),
           padj < padj_cutoff,
           log2FoldChange < -log2fc_cutoff) %>%
    pull(ensembl_id) %>%
    unique()

message(sprintf("Background genes: %d", length(background_genes)))
message(sprintf("Up-regulated genes: %d", length(up_genes)))
message(sprintf("Down-regulated genes: %d", length(down_genes)))

run_and_save_go <- function(genes, direction, output_csv, output_pdf) {
    if (length(genes) == 0) {
        message(sprintf("No significant %s genes. Skipping GO enrichment.", direction))
        write.csv(data.frame(), output_csv, row.names = FALSE)
        pdf(output_pdf)
        plot.new()
        text(0.5, 0.5, sprintf("No significant %s genes for GO", direction),
             cex = 1.2, col = "gray40")
        dev.off()
        return(invisible(NULL))
    }

    ego <- enrichGO(gene          = genes,
                    universe      = background_genes,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = "ENSEMBL",   
                    ont           = "ALL",        
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = FALSE)        

    ego_df <- as.data.frame(ego)
    write.csv(ego_df, output_csv, row.names = FALSE)

    if (nrow(ego_df) > 0) {
        title_text <- sprintf("GO: %s in %s KO vs WT", direction, comparison_name)
        p <- dotplot(ego, showCategory = 10, title = title_text)
        ggsave(filename = output_pdf, plot = p, device = "pdf",
               width = 8, height = 10, units = "in", dpi = 300)
    } else {
        message(sprintf("No significant GO terms for %s genes.", direction))
        pdf(output_pdf)
        plot.new()
        text(0.5, 0.5, "No significant GO terms found",
             cex = 1.2, col = "gray40")
        dev.off()
    }
}

run_and_save_go(up_genes,   "Upregulated",
                snakemake@output[["go_up"]],   snakemake@output[["plot_up"]])
run_and_save_go(down_genes, "Downregulated",
                snakemake@output[["go_down"]], snakemake@output[["plot_down"]])
