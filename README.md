# SNO_proteomics

**Differential Abundance of S-Nitrosylated Proteins and Functional Enrichment Analysis**

This repository contains the R code and data analysis pipeline used to assess differential abundance of S-nitrosylated (SNO) proteins in human liver cancer cells. The analysis compares parental (HepG2) and Sorafenib-resistant (R-HepG2) cell lines, under treated and untreated conditions, followed by functional enrichment using KEGG pathways.

## Overview

The main components of the analysis include:

* **Data normalization** using a class-specific quantile (CSQ) strategy to account for treatment and cell-type differences.
* **Differential abundance analysis** performed with the `limma` package on CSQ-normalized data.
* **Over-representation analysis (ORA)** of significantly changing proteins using KEGG pathway annotations.

## Repository Contents

* `SNO_proteomics.Rmd` — R Markdown file containing the full analysis pipeline.
* `SNO_proteomics.html` — Rendered output of the Rmd file.
* `SNO_proteomics.R` — R script used to run the pipeline and generate figures and tables separately.
* `functions_sno.R` — Custom helper functions used in the analysis.
* `SNO_proteomics.Rproj` — RStudio project file for convenient setup.
* `data/`` — Folder containing input data (proteomics results, translation table).
* `README.md` — This description of the project and usage.

## Reproducibility

This repository includes all code and data needed to reproduce the results shown in the manuscript. To replicate the analysis:

1. Open the RStudio project: `SNO_proteomics.Rproj`
2. Run or knit `SNO_proteomics.Rmd` to generate the full report.
3. Alternatively, source `SNO_proteomics.R` to generate figures and tables separately.

All required packages are installed/loaded in the Rmd file. A `sessionInfo()` is printed at the end for transparency on package versions.

