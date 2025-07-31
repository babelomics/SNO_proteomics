#!/usr/bin/env Rscript

# Purpose: DEA and KEGG pathway ORA of SNO-proteins
# Author: Miriam Payá-Milans
# Date: 31-07-2025

#########################
### Environment Setup ###
#########################

source("functions_sno.R")

# List required packages and install if missing from CRAN/Bioconductor

packages <- c("readxl", "openxlsx", "readr", "tibble", "dplyr", "stringr", "tidyr", 
              "limma", "clusterProfiler", "ggplot2", "ggvenn", "ggnewscale", "renv")
package_manager(packages, "install")

# Load required packages

suppressMessages(library(openxlsx))
suppressMessages(library(readr))
suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggvenn))
suppressMessages(library(limma))
suppressMessages(library(clusterProfiler))


####################
### Data Loading ###
####################

# Set working directory and Read data files. If already parsed, load and skip section
datadir <- "./data"
if (file.exists(file.path(datadir, "data_parsed.Rdata"))) {
  load(file.path(datadir, "data_parsed.Rdata"))
  message("Data loaded.")
} else {
  # Load data
  f <- "Proteomics_Results.xlsx"
  prot_res <- readxl::read_excel(file.path(datadir, f), .name_repair = "unique_quiet")
  tr_df <- readr::read_tsv(file.path(datadir, "UPacc2EG.tsv"), show_col_types = F)
  
  ## Extract Signal intensity as Relative Protein Abundance
  prot_abund <- prot_res %>%
    dplyr::select(Accession, contains("Area")) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("Accession") %>%
    dplyr::rename_with(\(x) stringr::str_extract(x, "(^.+):", group = 1))
  
  # remove proteins with >50% zero values (192)
  prot_abund <- prot_abund[rowSums(prot_abund == 0) < ncol(prot_abund)/2,]
  
  ## Extract Sample Metadata
  case_info <- Filter(\(x) stringr::str_detect(x, "Area"), colnames(prot_res))
  case_tr <- c(
    ".*2 control R.*" = "control_untr",
    ".*2 sorafenib R.*" = "control_sora",
    ".*2 resistentes R.*" = "resist_untr",
    ".*2 resistentes sorafenib R.*" = "resist_sora"
  )
  
  case_metadata <- tibble(Sample = stringr::str_extract(case_info, "(^.+):", group = 1),
                          Condition = stringr::str_replace_all(case_info, case_tr)) %>%
    tidyr::separate_wider_delim(Condition, "_", cols_remove = F,
                                names = c("Cell_line", "Treated")) %>%
    dplyr::mutate(Condition = factor(Condition, levels = case_tr),
                  Cell_line = factor(Cell_line, levels = c("control", "resist")),
                  Treated = factor(Treated, levels = c("untr", "sora")))
  
  ## Extract Protein Metadata and Rename columns
  protein_metadata <- prot_res %>%
    dplyr::mutate(
      Uniprot_ID = stringr::str_extract(.$Description, "PE.+\\[(.*?)\\]", group = 1),
      Protein_name = stringr::str_extract(.$Description, "(^.+) OS=", group = 1)
    ) %>%
    left_join(dplyr::select(tr_df, from, symbol, gene_id), by = c("Accession" = "from")) %>%
    dplyr::select(
      Protein_name,
      Accession,
      Uniprot_ID,
      Gene_id = gene_id,
      Symbol = symbol,
      Coverage = `ΣCoverage`,
      n_Uniq_Pept = `Σ# Unique Peptides`,
      n_AAs = `# AAs`,
      MW_kDa = `MW [kDa]`,
      calc_pI = `calc. pI`
    ) %>%
    dplyr::mutate(MW_kDa = round(MW_kDa, 3), calc_pI =  round(calc_pI, 2)) %>%
    dplyr::distinct() %>% 
    mutate(Protein_name = str_remove(Protein_name, "^.* similar to "), 
           Protein_name = str_remove(Protein_name, ", mRNA$"), 
           Protein_name = str_remove(Protein_name, "Homo sapiens "), 
           Protein_name = str_remove(Protein_name, " \\(Fragment\\)"))
  
  save(prot_abund, case_metadata, protein_metadata, 
       file = file.path(datadir, "data_parsed.Rdata"))
}


#####################
### Normalization ###
#####################
message("Running Normalization.")

# create dir for results if needed
normdir <- "./01_normalization"
if (!dir.exists(normdir)){dir.create(normdir)}

## Data Normalization
# 1. Replace zeros with NA and Apply log2 transformation
log2_data <- prot_abund
log2_data[log2_data == 0] <- NA
log2_data <- log2(log2_data)

# 2. Class-specific quantile normalization
conditions <- unique(case_metadata$Condition)
norm_list <- list()
for (cond in conditions) {
  cond_samples <- case_metadata %>% filter(Condition == cond) %>% pull(Sample)
  sub_data <- log2_data[, cond_samples, drop = FALSE]
  norm_data <- normalizeBetweenArrays(sub_data, method = "quantile")
  norm_list[[cond]] <- norm_data
}
norm_class_quantile <- do.call(cbind, norm_list) %>% as.data.frame()

write_tsv(norm_class_quantile %>% rownames_to_column("Accession"), 
          file.path(normdir, "normalized_data.tsv"), na = "")
# table(rowSums(is.na(norm_class_quantile)))

## VISUALIZATION ##
## Wrangling of Normalized values for visualization
# Combine all datasets into a long-format dataframe
preplot_data <- bind_rows(
  log2_data %>%
    tibble::rownames_to_column("Accession") %>%
    tidyr::pivot_longer(-Accession, names_to = "Sample", values_to = "Area") %>%
    dplyr::mutate(Method = "Raw (log2)"),
  
  norm_class_quantile %>%
    tibble::rownames_to_column("Accession") %>%
    tidyr::pivot_longer(-Accession, names_to = "Sample", values_to = "Area") %>%
    dplyr::mutate(Method = "Class-sp Quantile")
)

# clean NA values and Add sample condition metadata
plot_data <- preplot_data %>%
  dplyr::filter(!is.na(Area)) %>%
  dplyr::mutate(Method = factor(Method,
                                levels = unique(preplot_data$Method))) %>%
  dplyr::left_join(case_metadata, by = "Sample") %>%
  dplyr::left_join(protein_metadata, by = "Accession") %>%
  dplyr::mutate(avg = mean(Area), .by = c("Method", "Condition"))

## Combined visualization of normalization results
# Violin plots of samples faceted by method
outfile <- file.path(normdir, "norm_boxplot.png")
p <- plot_violin_box(plot_data, xvar = Sample, yvar = Area, fill_var = Condition, 
                     facet_var = "Method", axis_scale = "free_y", 
                     title = "Normalization comparison by Sample", 
                     file_path = outfile, w = 4, h = 6)

# Density plots of Samples faceted by method
outfile <- file.path(normdir, "norm_density.png")
p <- plot_density(plot_data, Area, group_var = Sample, color_var = Condition, 
                  facet_var = "Method", axis_scale = "free_x", 
                  title = "Protein abundance distribution:\n Samples colored by Condition", 
                  file_path = outfile, w = 4, h = 6)

# Density plots of Condition-aggregated samples faceted by method
outfile <- file.path(normdir, "norm_density_grouped.png")
p <- plot_density(plot_data, Area, color_var = Condition, 
                  facet_var = "Method", axis_scale = "free_x", 
                  title = "Protein abundance distribution:\n Samples averaged by Condition", 
                  file_path = outfile, w = 4, h = 6)

# Violin plots of Proteins faceted by method
n_prot <- 25
p <- plot_violin_box(plot_data, xvar = Accession, yvar = Area, 
                     max_terms = n_prot, facet_var = "Method", 
                     title = "Intensity of Protein Accessions across samples")
p <- p + geom_point(aes(y = avg, color = Condition), shape = 23, stroke = 1, size = .5) + 
  labs(color = "Mean Intensity\nby Condition")
ggsave(file.path(normdir, paste0("protein_boxplot_", n_prot, ".png")), p, 
       width = 6, height = 8)

# Density plots of Proteins faceted by method
n_prot <- 20
outfile <- file.path(normdir, paste0("protein_density_", n_prot, ".png"))
p <- plot_density(plot_data, Area, group_var = Accession, max_group = n_prot,  
                  facet_var = "Method", lw = 0.3, axis_scale = "fixed", 
                  title = "Distribution of Protein Intensities across samples", 
                  file_path = outfile, w = 5, h = 7)

# PCA of all methods
outfile <- file.path(normdir, "norm_pca_plot.png")
p <- run_plot_pca(plot_data, run_by = "Method", plot_group = Condition,
                  title = "PCA of SNO dataset; raw (log2) and normalized data.", 
                  file_path = outfile, w = 5.5, h = 4)


##################################
### Testing of group Contrasts ###
##################################
message("Running Testing of group Contrasts.")

# Create output folder if needed
diffdir <- "02_diff_analysis"
if (!dir.exists(diffdir)) dir.create(diffdir, showWarnings = FALSE)

## RUN DIFFERENTIALLY STATISTICAL ANALYSIS ##
# Run automatically-designed contrasts on limma; add protein annotations
limma_results <- fit_contrasts(norm_class_quantile, case_metadata$Condition) %>% 
  dplyr::left_join(protein_metadata, by = "Accession")

## TABLES ##
# number of significant proteins
pval <- 0.05
limma_results %>% 
  summarize(n_pval = sum(P.Value < pval, na.rm = T), 
            n_adj = sum(adj.P.Val < 0.05, na.rm = T), .by = "Contrast") %>% 
  write_tsv(file.path(diffdir, "n_signif_res.tsv"))

# soft filtering by p.value < 0.05, in view of poor results with padj
limma_results_sig <- limma_results %>%
  dplyr::filter(P.Value < pval)

# Write limma results to multiple sheets
save_xlsx(limma_results, "Contrast", 
          file.path(diffdir, "limma_results_full.xlsx") )
save_xlsx(limma_results_sig, "Contrast", 
          file.path(diffdir, sprintf("limma_results_p%s.xlsx", pval)) )


## VISUALIZATION ##
# Calculation of combined score to include protein detection as coverage x mass
limma_results_sig <- limma_results_sig %>% 
  dplyr::mutate(n_contrast = length(Contrast), .by = Accession) %>% 
  dplyr::mutate(abs_Change = abs(t), 
                log10_pvalue = -log10(P.Value), 
                CombinedScore = abs(t)*log10_pvalue*n_contrast*Coverage*MW_kDa)

# Balloon plot of Proteins vs Contrasts
max_terms <- 20
p <- plot_balloon(limma_results_sig, xvar = Contrast, yvar = Accession, 
                  max_y = max_terms, color_var = logFC, size_var = log10_pvalue, 
                  size_lab = "-log10(p)", sort_var = "CombinedScore", outdir = diffdir, 
                  prefix = "limma", p_type = "P.Value", pval = pval, ylab = "Protein", 
                  title = "Differential Abundance Balloon Plot (p<0.05)")

# alter the balloon plot for limma results to add annotations
p <- p + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

# select row annotation data to plot
annotation_data <- protein_metadata %>%
  dplyr::select(Accession, calc_pI, MW_kDa) %>%
  dplyr::filter(Accession %in% unique(p$data$Accession)) %>% 
  dplyr::mutate(Accession = factor(Accession, levels = rev(unique(p$data$Accession)))) %>% 
  tidyr::pivot_longer(cols = c(calc_pI, MW_kDa), names_to = "Feature", values_to = "Value")

m1 <- annotation_data %>% 
  dplyr::filter(Feature == "calc_pI")
m2 <- annotation_data %>% 
  dplyr::filter(Feature == "MW_kDa")

# plot annotation columns with 2 color scales
p_annot <-ggplot(mapping = aes(x = Feature, y = Accession)) +
  # This bit is for making scales
  geom_tile(data=m1, aes(fill = Value)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 7) +
  guides(fill = guide_legend(title="calc_pI")) +
  ggnewscale::new_scale_fill() +
  geom_tile(data=m2, aes(fill=Value)) +
  scale_fill_gradient(low = "white", high = "black") + 
  guides(fill = guide_legend(title="MW_kDa")) +
  theme_minimal() +
  labs(y = "Protein") + 
  scale_x_discrete(position = "top") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        axis.title.x = element_blank(),
        axis.text.y = element_text(),
        legend.position = "top", 
        legend.box = "vertical", 
        legend.key.spacing = unit(0.1, 'cm'),
        legend.spacing = unit(0, 'cm'), 
        legend.key.size = unit(0.5, 'cm'))

# combine plots
comb_plot <- p_annot + p + patchwork::plot_layout(widths = c(0.2, 1))
outfile <- sprintf("balloon_limma-n%d-p%g-row_annot.png", max_terms, pval)
ggsave(file.path(diffdir, outfile), comb_plot, 
       width = 7, height = 7)

# extract top protein names
top_ids <- p$data %>% 
  dplyr::select(Accession, Protein_name, Gene_id, Symbol) %>% distinct()
outfile <- sprintf("balloon_limma-n%d-p%g-top_ids.tsv", max_terms, pval)
write_tsv(top_ids, file.path(diffdir, outfile))


#############################
### Functional Enrichment ###
#############################
## Functional Enrichment Analysis with CLUSTERPROFILER
message("Running Functional Enrichment Analysis.")

# Set and Create output folder if missing
kegg_dir <- "03_func_annot"
if (!dir.exists(kegg_dir)) dir.create(kegg_dir, showWarnings = FALSE)

# Extract significant proteins per contrast
contrast_list <- unique(limma_results_sig$Contrast)
sig_prot <- list()
for (c in contrast_list) {
  sig_prot[[c]] <- limma_results_sig %>% dplyr::filter(Contrast == c) %>% pull(Accession)
}

# Translate into significant genes and apply KEGG GSEA
kk <- list()
for (contrast_name in contrast_list) {
  # Extract significant genes for this contrast
  contrast_genes <- protein_metadata %>%
    dplyr::filter(Accession %in% sig_prot[[contrast_name]]) %>%
    dplyr::pull(Gene_id) %>% na.omit() %>% unique()
  
  # Run enrichment with whole universe
  kk[[contrast_name]] <- enrichKEGG(gene = contrast_genes,
                                        organism = "hsa",
                                        pvalueCutoff = 0.05)
}

# Extract and bind the resulting outputs
keggres <- lapply(seq_along(contrast_list), function(i) {
  kk[[i]]@result %>%
    dplyr::mutate(Contrast = contrast_list[i], .before = 1)
}) %>% 
  dplyr::bind_rows(.) %>%
  dplyr::mutate(symbolID = purrr::map_chr(geneID, mapIDs, map_df = protein_metadata),
                Accession = purrr::map_chr(geneID, mapAcc, map_df = protein_metadata)) %>%
  dplyr::mutate(n_contrast = sum(pvalue < 0.05), .by = Description) %>% 
  dplyr::mutate(log10_pvalue = round(-log10(pvalue), 3),
                CombinedScore = FoldEnrichment * log10_pvalue * zScore * n_contrast)

## Write KEGG results to multiple sheets into one file
save_xlsx(keggres, "Contrast", file.path(kegg_dir, "kegg_results.xlsx"))

## VISUALIZATION ##
# Balloon plot: by p-adjust < 0.05
max_terms <- Inf
con_prefix <- "KEGG_contrasts"
p <- plot_balloon(keggres, xvar = Contrast, yvar = Description, 
                  max_y = max_terms, size_var = Count, color_var = FoldEnrichment,
                  p_type = "p.adjust", sort_var = "CombinedScore", 
                  ylab = "KEGG Pathway", prefix = con_prefix, outdir = kegg_dir, 
                  title = "KEGG Enrichment on significant genes (p.adjust)")

# Balloon plot: kegg vs prot, by p-adjust < 0.05
keggres_gene <- keggres %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  separate_longer_delim(symbolID, "/") %>% 
  dplyr::filter(symbolID != "NA")
p <- plot_balloon(keggres_gene, xvar = symbolID, yvar = Description, 
                  max_y = max_terms, size_var = Count, color_var = FoldEnrichment, 
                  p_type = "p.adjust", xlab = "Gene Symbol", ylab = "KEGG Pathway", 
                  prefix = con_prefix, outdir = kegg_dir, sort_var = "CombinedScore", 
                  title = "KEGG Enrichment across contrasts with participating genes")

# extract plotted protein names
pids <- p$data %>% 
  dplyr::select(symbolID) %>% distinct() %>% 
  left_join(protein_metadata %>% 
              dplyr::select(Accession, Protein_name, Gene_id, Symbol),
            by = c("symbolID" = "Symbol"))
outfile <- sprintf("balloon_%s-%dgenes-p%g.tsv", con_prefix, length(unique(pids$symbolID)), pval)
write_tsv(pids, file.path(kegg_dir, outfile))

# Balloon plot: by p-value < 0.05
p <- plot_balloon(keggres, xvar = Contrast, yvar = Description, 
                  max_y = max_terms, size_var = Count, color_var = FoldEnrichment,
                  p_type = "pvalue", sort_var = "CombinedScore", 
                  ylab = "KEGG Pathway", prefix = con_prefix, outdir = kegg_dir, 
                  title = "KEGG Enrichment on significant genes (pvalue)")

# Balloon plot: by p-value < 0.05, removing disease terms.
p <- plot_balloon(keggres %>% dplyr::filter(!str_detect(category, "Disease")), 
                  xvar = Contrast, yvar = Description, 
                  max_y = max_terms, size_var = Count, color_var = FoldEnrichment,
                  p_type = "pvalue", sort_var = "CombinedScore", 
                  ylab = "KEGG Pathway", prefix = "KEGG_noDisease", outdir = kegg_dir, 
                  title = "KEGG Enrichment on significant genes (pvalue)")

# Dotplots from Iteration over KEGG results in the list
plot_dot(keggres, "Description", "p.adjust", prefix = "KEGG", outdir = file.path(kegg_dir, "dotplots"))
plot_dot(keggres, "Description", "pvalue", prefix = "KEGG", outdir = file.path(kegg_dir, "dotplots"))


# specific plots by contrast
# Balloon plot: C_sora vs C_untr, kegg vs prot, by pvalue < 0.05
csvscu <- "control_sora-control_untr"
keggres_csvscu <- keggres %>% 
  dplyr::filter(Contrast == csvscu, pvalue < 0.05) %>% 
  separate_longer_delim(cols = c(geneID, symbolID, Accession), delim = "/") %>% 
  dplyr::filter(symbolID != "NA") %>% 
  inner_join(dplyr::filter(limma_results_sig, Contrast == csvscu)[,1:4], by = "Accession") %>% 
  mutate(symbolID = factor(symbolID, levels = unique(.$symbolID)), 
         FCxFE = logFC * FoldEnrichment, 
         FCxFE2 = abs(logFC) * FoldEnrichment, 
         FCxFE3 = -logFC * FoldEnrichment, 
         sortScore = FoldEnrichment * log10_pvalue * zScore)
p2 <- plot_balloon(keggres_csvscu, xvar = symbolID, yvar = Description, 
                   max_y = 30, size_var = FCxFE, color_var = logFC, 
                   sort_var = "sortScore", 
                   prefix = "KEGG_CSvsCU", outdir = kegg_dir, 
                   p_type = "pvalue", xlab = "Gene Symbol", ylab = "KEGG Pathway", 
                   title = "KEGG Enrichment on control cells")

# extract plotted protein names
p2$data %>% 
  dplyr::select(symbolID) %>% distinct() %>% 
  left_join(protein_metadata %>% 
              dplyr::select(Accession, Protein_name, Gene_id, Symbol),
            by = c("symbolID" = "Symbol")) %:% 
  write_tsv(file.path(kegg_dir, "balloon_KEGG_CSvsCU-prot.tsv"))

# Balloon plot: R_sora vs C_sora, kegg vs prot, by pvalue < 0.05
rsvscs <- "resist_sora-control_sora"
keggres_rsvscs <- keggres %>% 
  dplyr::filter(Contrast == rsvscs, pvalue < 0.05) %>% 
  separate_longer_delim(cols = c(geneID, symbolID, Accession), delim = "/") %>% 
  dplyr::filter(symbolID != "NA") %>% 
  inner_join(dplyr::filter(limma_results_sig, Contrast == rsvscs)[,1:4], by = "Accession") %>% 
  mutate(symbolID = factor(symbolID, levels = unique(.$symbolID)), 
         FCxFE = logFC * FoldEnrichment, 
         FCxFE2 = abs(logFC) * FoldEnrichment, 
         FCxFE3 = -logFC * FoldEnrichment, 
         sortScore = FoldEnrichment * log10_pvalue * zScore)
p3 <- plot_balloon(keggres_rsvscs, xvar = symbolID, yvar = Description, 
                   max_y = 30, size_var = FCxFE3, size_lab = "-FCxFE", 
                   color_var = logFC, sort_var = "sortScore", 
                   prefix = "KEGG_RSvsCS", outdir = kegg_dir, 
                   p_type = "pvalue", xlab = "Gene Symbol", ylab = "KEGG Pathway", 
                   title = "KEGG Enrichment on Sorafenib-treated cells")

# extract plotted protein names
p3$data %>% 
  dplyr::select(symbolID) %>% distinct() %>% 
  left_join(protein_metadata %>% 
              dplyr::select(Accession, Protein_name, Gene_id, Symbol),
            by = c("symbolID" = "Symbol")) %:% 
  write_tsv(file.path(kegg_dir, "balloon_KEGG_RSvsCS-prot.tsv"))


####################
### Session Info ###
####################
message("Writing Session Info")

writeLines(capture.output(sessionInfo()), "session_info.txt")
