#!/usr/bin/env Rscript


# Functions to install/update/remove packages from CRAN/Bioconductor

package_manager <- function(pkgs, action = c("install", "remove", "update"), ask = FALSE) {
  action <- match.arg(action)
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  switch(action,
         install = {
           for (pkg in pkgs) {
             if (!requireNamespace(pkg, quietly = TRUE)) {
               if (pkg %in% BiocManager::available()) {
                 BiocManager::install(pkg, ask = ask)
               } else {
                 install.packages(pkg)
               }
               message("Installed: ", pkg)
             } else {
               # message("Already installed: ", pkg)
             }
           }
         },
         remove = {
           renv::remove(pkgs)
           renv::clean()
           message("Removed: ", paste(pkgs, collapse = ", "))
         },
         update = {
           for (pkg in pkgs) {
             if (ask) {
               ans <- readline(paste("Update", pkg, "? [y/N]: "))
               if (tolower(ans) != "y") next
             }
             if (pkg %in% BiocManager::available()) {
               BiocManager::install(pkg, update = TRUE, ask = ask)
             } else {
               install.packages(pkg)
             }
             message("Updated: ", pkg)
           }
         }
  )
  renv::snapshot()
  renv::status()
}


# Function to generate violin + boxplot + avg points
plot_violin_box <- function(data, xvar, yvar, fill_var = NULL, max_terms = Inf, 
                            facet_var = ".", axis_scale = "fixed", ncol = 1,
                            title, file_path = NULL, w = 6, h = 6) {
  x_chr <- deparse(substitute(xvar))
  if (max_terms < length(unique(data[[x_chr]]))) {
    data <- data %>% filter(.data[[x_chr]] %in% unique(.data[[x_chr]])[1:max_terms])
  }
  p <- ggplot(data, aes(x = {{xvar}}, y = {{yvar}}, fill = {{fill_var}})) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2) +
    geom_point(aes(y = .data$avg), shape = 4, size = 1, color = "black") +
    facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol, scales = axis_scale) +
    # facet_grid(as.formula(paste(facet_var, "~ .")), scales = axis_scale) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.x = if (facet_var == "."){element_blank()} else {element_text()}, 
          legend.position = "top") + 
    guides(fill = guide_legend(nrow = 2,byrow = TRUE)) +
    scale_fill_manual(values = c('red', 'green', 'blue', "orange")) +
    labs(title = title, y = "Log2 Intensity", x = x_chr)
  if (!is.null(file_path)){
    ggsave(file_path, p, width = w, height = h)
  }
  return(p)
}

# Function to generate density plots
plot_density <- function(data, xvar, group_var = NULL, color_var = NULL, max_group = Inf, 
                         facet_var = ".", ncol = 1, axis_scale = "fixed", lw = 0.5,
                         title, file_path = NULL, w = 6, h = 6) {
  group_chr <- deparse(substitute(group_var))
  if(!is.null(group_chr)) 
    if (max_group < length(unique(data[[group_chr]]))) {
      data <- data %>% filter(.data[[group_chr]] %in% unique(.data[[group_chr]])[1:max_group])
    }
  p <- ggplot(data, aes(x = {{xvar}}, group = {{ group_var }}, color = {{color_var}})) +
    geom_density(aes(linetype = {{ color_var }}), linewidth = lw) +
    facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol, scales = axis_scale) +
    theme_bw() +
    theme(strip.text.x = if (facet_var == ".") element_blank() else element_text(), 
          legend.position = "top") +
    guides(color = guide_legend(nrow = 2,byrow = TRUE)) +
    scale_color_manual(values = c('red', 'green', 'blue', "orange")) +
    labs(title = title, x = "Log2 Intensity", y = "Density")
  if (!is.null(file_path)){
    ggsave(file_path, p, width = w, height = h)
  }
  return(p)
}

# Auxiliary function to run PCA on a dataset
# NA values are converted to 0
run_pca <- function(df, method_label) {
  # Wide format: proteins as rows, samples as columns
  mat <- df %>%
    filter(Method == method_label) %>%
    select(Accession, Sample, Area) %>%
    pivot_wider(names_from = Sample, values_from = Area) %>%
    column_to_rownames("Accession") %>%
    as.matrix()
  
  mat <- t(mat)  # PCA on samples
  mat[is.na(mat)] <- 0  # remove NA values
  
  pca <- prcomp(mat, scale. = TRUE)
  # % variance explained
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
  
  pca_df <- as.data.frame(pca$x[, 1:2])  # First 2 PCs
  pca_df$Sample <- rownames(pca_df)
  pca_df$Method <- method_label
  pca_df$Condition <- df$Condition[match(pca_df$Sample, df$Sample)]
  pca_df$FacetLabel <- factor(paste0(
    method_label, "\nPC1: ", var_exp[1], "%, PC2: ", var_exp[2], "%"
  ))
  return(pca_df)
}

# Function to call the full PCA pipeline
run_plot_pca <- function(df, run_by = "Method", plot_group = NULL,
                         title, file_path = NULL, w = 6, h = 8){
  # group_chr <- deparse(substitute(plot_group))
  methods <- unique(df[[run_by]])
  
  pca_res <- bind_rows(lapply(methods, function(m) run_pca(df, m)))
  
  p <- ggplot(pca_res, aes(x = PC1, y = PC2, fill = {{plot_group}}, shape = {{plot_group}})) +
    geom_point(size = 3) +
    facet_wrap(~ FacetLabel, scales = "free") +
    labs(title = title) + 
    scale_fill_manual(values = c('red', 'green', 'blue', "orange")) +
    scale_shape_manual(values = c(22,23,24,25)) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold", size = 12),
          legend.position = "bottom")
  if (!is.null(file_path)){
    ggsave(file_path, p, width = w, height = h)
  }
  return(p)
}



# Function for Differential Protein Abundance Analysis with Limma
fit_contrasts <- function(norm_data, case_factor = F, desmat = F, conmat = F) {
  if (!isFALSE(case_factor)) {
    # 1. Define the design matrix
    design <- model.matrix(~0 + case_factor)
    colnames(design) <- levels(case_factor)
    
    # 2. Create all pairwise contrasts
    conds <- levels(case_factor)
    contrast_list <- combn(rev(conds), 2, FUN = function(x) paste0(x[1], "-", x[2]))
    contrast.matrix <- makeContrasts(contrasts = contrast_list, levels = design)
  } else if (!isFALSE(desmat) & !isFALSE(conmat)){
    design <- desmat
    contrast.matrix <- conmat
  } else {warning("Either case factor or desmat+conmat missing.")}
  
  # 3. Fit the model and apply contrasts
  fit <- lmFit(norm_data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # 4. Extract all top tables and create a single dataframe with all contrasts
  tt_list <- lapply(seq_along(contrast_list), function(i) {
    topTable(fit2, coef = i, number = nrow(norm_data), adjust.method = "BH", sort.by = "P") %>%
      tibble::rownames_to_column("Accession") %>%
      mutate(Contrast = contrast_list[i])
  })
  all_results <- dplyr::bind_rows(tt_list)
  
  return(all_results)
}

# Function to map gene IDs to symbols
mapIDs <- function(gene_ids_string, map_df) {
  ids <- unlist(str_split(gene_ids_string, "/"))
  symbols <- map_df$Symbol[match(ids, map_df$Gene_id)]
  paste(symbols, collapse = "/")
}

# Function to map gene IDs to protein accessions
mapAcc <- function(gene_ids_string, map_df) {
  ids <- unlist(str_split(gene_ids_string, "/"))
  prot_acc <- map_df$Accession[match(ids, map_df$Gene_id)]
  paste(prot_acc, collapse = "/")
}

## Visualization of significant results; filtering by CombinedScore
# BALLOON PLOT

get_top_ids <- function(diffres, yvar, max_y = Inf, xvar = NULL, max_x = Inf, 
                        sort_var = "CombinedScore", p_type = "p.adjust", pval = 0.05){
  top_diff <- diffres %>%
    dplyr::filter(.data[[p_type]] < pval) %>%
    group_by(.data[[yvar]]) %>%
    dplyr::mutate(mean_sortvar = mean(.data[[sort_var]])) %>%
    ungroup() %>%
    dplyr::arrange(desc(mean_sortvar)) %>%
    dplyr::filter(.data[[yvar]] %in% unique(.data[[yvar]])[1:max_y])
  if (!is.null(xvar)) {
    top_diff <- top_diff %>%
      dplyr::filter(.data[[xvar]] %in% unique(.data[[xvar]])[1:max_x])
  }
  top_diff <- top_diff %>%
    dplyr::mutate(colvar = factor(.data[[yvar]], levels = rev(unique(.data[[yvar]]))))
  return(top_diff)
}


plot_balloon <- function(diffres, xvar, yvar, color_var = NULL, shape_var = NULL, 
                         size_var = NULL, size_lab = NULL, 
                         p_type = "p.adjust", pval = 0.05, 
                         max_x = Inf, max_y = Inf, sort_var = p_type, 
                         title = "Balloon plot", xlab = NULL, ylab = NULL, 
                         outdir = NULL, prefix = "plot"){
  # extract string names from variables
  x_str <- deparse(substitute(xvar))
  y_str <- deparse(substitute(yvar))
  xlab <- if (is.null(xlab)) x_str else xlab
  ylab <- if (is.null(ylab)) y_str else ylab
  max_x <- min(max_x, length(unique(diffres[[x_str]])))
  max_y <- min(max_y, length(unique(diffres[[y_str]])))
  size_lab <- if (is.null(size_lab)) deparse(substitute(size_var)) else size_lab
  
  # crop/filter/sort the data to max items per axis and significance
  top_diff <- get_top_ids(diffres, y_str, max_y, x_str, max_x, sort_var, p_type, pval)
  if (nrow(top_diff) == 0) stop("No significant results.")
  
  # prepare plot
  p <- ggplot(top_diff, aes(x = {{xvar}}, y = colvar)) + 
    geom_point(aes(color = {{color_var}}, shape = {{shape_var}}, 
                   size = {{size_var}})) + 
    scale_size_continuous(name = size_lab) + 
    scale_shape_manual(values = seq(15,18)) + 
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_x_discrete(position = "top") +
    labs(title = title, x = xlab, y = ylab) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0),
          plot.title = element_text(hjust = 0.6), 
          legend.spacing = unit(0, 'cm'), 
          legend.key.size = unit(0.5, 'cm'))

  # Save plot
  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)
    n_cols <- length(unique(p$data[[x_str]]))
    n_rows <- length(unique(p$data[[y_str]]))
    max_x_char <- max(nchar(as.character(p$data[[x_str]])))
    max_y_char <- max(nchar(as.character(p$data[[y_str]])))
    max_labs <- max(c(nchar(xlab), nchar(ylab), nchar(size_lab)))
    w <- 1.5 + max_y_char * 0.05 + n_cols * 0.2 + max_labs * 0.1
    h <- 1 + max_x_char * 0.05 + n_rows * 0.2
    outfile <- file.path(outdir, sprintf("balloon_%s-%s_vs_%s-%s%g.png", 
                                         prefix, xlab, ylab, p_type, pval))
    ggsave(outfile, plot = p, width = w, height = h)
    
    # Save KEGG data
    if ("FoldEnrichment" %in% colnames(p$data)) { 
      p_data <- p$data %>% 
        dplyr::select(Contrast, category, Description, FoldEnrichment, pvalue, 
                      p.adjust, geneID, symbolID, Accession) %>% 
        distinct()
      write_tsv(p_data, stringr::str_replace(outfile, "png", "tsv"))
      
      path_data <- distinct(dplyr::select(p$data, category, subcategory, Description))
      write_tsv(path_data, stringr::str_replace(outfile, ".png", "-path.tsv"))
    }
  }
  
  return(p)
}


# DOTPLOT: Iterate over KEGG results in the list
plot_dot <- function(keggres, yvar, p_type = "p.adjust", prefix = NULL,
                     max_terms = 25, pval = 0.05, outdir){
  # try to create outdir on existing folders
  if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)

  plot_by <- "Contrast"
  for (name in unique(keggres[[plot_by]])) {
    enrich_df <- keggres %>%
      dplyr::filter(.data[[plot_by]] == name)
    aux <- get_top_ids(enrich_df, yvar, max_terms, p_type = p_type, pval = pval)
    if (nrow(aux) == 0) next
    
    # Create the plot
    p <- ggplot(aux, aes(x = CombinedScore, y = colvar, size = Count, fill = log10_pvalue)) +
      geom_point(shape = 21) +
      labs(title = sprintf("KEGG enrichment (%s)", p_type),
           subtitle = name, x = "Combined Score", y = "KEGG Pathway") +
      scale_fill_gradient(low = "khaki1", high = "firebrick2",
                          labels = scales::number_format(accuracy = 0.1)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
      guides(size = "none") +
      theme_light() +
      theme(plot.title = element_text(hjust = 1, size = 12),
            plot.subtitle = element_text(hjust = 1, size = 10))
    
    # Save plot and table
    w <- 1000 + max(nchar(as.character(aux$Description))) * 15
    h <- 1000 + nrow(aux) * 15
    outfile <- paste0(paste(c("dotplot",prefix, p_type, name), collapse = "_"), ".png")
    ggsave(file.path(outdir, outfile), plot = p, width = w, height = h, units = "px")
    # readr::write_tsv(aux, file.path(outdir, sprintf("dotplot_%s_%s_%s.tsv",prefix, p_type, name)))
  }
}

# Save to Excel
save_xlsx <- function(res, to_sheet, outname){
  wb <- createWorkbook()
  for (sheet in unique(res[[to_sheet]])) {
    addWorksheet(wb, sheet)
    res_filt <- dplyr::filter(res, .data[[to_sheet]] == sheet)
    writeData(wb, sheet, res_filt)
  }
  saveWorkbook(wb, outname, overwrite = TRUE)
}
