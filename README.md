# MethScope

<p align="center">
  <img src="docs/logo.png" alt="MethScope logo" width="200">
</p>

**MethScope** is an R package for ultra-fast analysis of sparse DNA methylome data using **Most Recurrent Methylation Patterns (MRMPs)**.

It supports downstream analysis for:
- Cell type annotation
- Cell type deconvolution
- Unsupervised clustering
- Cancer cell-of-origin prediction
- Missing value imputation

## Why MethScope?

Sparse single-cell and spatial methylome data are often too sparse to analyze directly. MethScope compresses methylation signals into MRMP-based embeddings so you can run robust and scalable downstream tasks with standard analysis workflows.

## Method overview

<p align="center">
  <img src="man/figures/overview.png" alt="MethScope workflow overview" width="980">
</p>

MethScope converts high-dimensional methylation atlas signals into compact MRMP features and applies these features across multiple analysis tasks.

Core workflow:
- Binarize methylation atlas profiles and consolidate recurrent patterns
- Select top recurrent methylation patterns (MRMPs)
- Encode each sample, cell, or pixel into an MRMP-based representation
- Run downstream modeling for annotation, deconvolution, imputation, and representation learning

Use cases supported in the current pipeline:
- Cell-type annotation in sparse single-cell methylome profiles
- Mini-bulk deconvolution for mixed-cell samples
- Missing-value imputation for sparse CpG measurements
- Representation learning for clustering and embedding analysis

## System Requirements

### Software Dependencies

- **R** >= 4.0 (tested on R 4.5.2 and R 4.6.0)
- **Operating systems tested**: macOS, Linux (Ubuntu)
- R package dependencies (installed automatically): `xgboost`, `dplyr`, `tidyr`, `stringr`, `caret`, `doParallel`, `parallel`, `ggplot2`, `uwot`, `magrittr`, `FNN`, `data.table`, `nnls`

### Hardware Requirements

No non-standard hardware is required. MethScope runs on a standard laptop or desktop CPU. No GPU is needed.

## Installation

Install from CRAN:

```r
install.packages("MethScope")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("zhou-lab/MethScope")
```

**Typical install time**: approximately 1 minute on a standard laptop.

## Demo

A small example dataset is bundled with the package in `inst/extdata/`. The following demo runs end-to-end cell type annotation using the included example `.cg` file and a pre-built mouse brain MRMP reference.

```r
library(MethScope)

# Locate bundled example files
example_file      <- system.file("extdata", "example.cg", package = "MethScope")
reference_pattern <- system.file("extdata", "Liu2021_MouseBrain.cm", package = "MethScope")

# Step 1: Generate cell-by-MRMP embedding matrix
input_pattern <- GenerateInput(example_file, reference_pattern)

# Step 2: Predict cell types using the built-in pre-trained mouse brain model
prediction_result <- PredictCellType(MethScope:::Liu2021_MouseBrain_P1000, input_pattern)

# Step 3: Visualize results
umap_plot <- PlotUMAP(input_pattern, prediction_result)
```

Expected output: a cell-by-MRMP matrix (`input_pattern`) and a data frame of predicted cell type labels with confidence scores (`prediction_result`). The UMAP plot will display cells colored by predicted cell type. Expected runtime on the bundled example: a few seconds.

## Tutorials and documentation

- Documentation website: [zhou-lab.github.io/MethScope](https://zhou-lab.github.io/MethScope/)
- End-to-end tutorial: [MethScope-Tutorial](https://zhou-lab.github.io/MethScope/articles/MethScope-Tutorial.html)
- Building MRMP references: [MethScope-MRMP](https://zhou-lab.github.io/MethScope/articles/MethScope-MRMP.html)

## Agent skill

This repository includes a reusable MethScope agent skill under `agent-skills/methscope/`.

### Codex

Install the skill under `$CODEX_HOME/skills/methscope/` with this layout:

```text
$CODEX_HOME/skills/methscope/
  SKILL.md        <- copy from agent-skills/methscope/codex/SKILL.md
  core/
```

Copy:
- `agent-skills/methscope/codex/SKILL.md` to `$CODEX_HOME/skills/methscope/SKILL.md`
- `agent-skills/methscope/core/` to `$CODEX_HOME/skills/methscope/core/`

Then invoke the skill when working on MethScope package usage, vignettes, `.cg` and `.cm` inputs, MRMP embeddings, prediction, training, deconvolution, or visualization.

### Claude

If you already keep a repository-level `CLAUDE.md`, copy the contents of `agent-skills/methscope/claude/CLAUDE.md` into it or reference that file from your existing Claude project instructions.

If you do not already have a project `CLAUDE.md`, use `agent-skills/methscope/claude/CLAUDE.md` as the starting project context for this repository.

Keep `agent-skills/methscope/core/` in the repository, because the Claude instructions point to those shared files.

### Shared files

- `agent-skills/methscope/core/INSTRUCTIONS.md`
- `agent-skills/methscope/core/WORKFLOWS.md`
- `agent-skills/methscope/core/REFERENCES.md`

## Data resources

- Example and reference data: [zhou-lab/methscope_data](https://github.com/zhou-lab/methscope_data)
- `.cg` generation and preprocessing: [YAME](https://zhou-lab.github.io/YAME/)
- Pattern interpretation: [knowYourCG](https://www.bioconductor.org/packages/release/bioc/html/knowYourCG.html)

## License

This project is licensed under the [GNU Affero General Public License v3.0 (AGPL-3.0)](LICENSE.md).

Copyright (c) 2025 Hongxiang Fu and Wanding Zhou (zhouw3@chop.edu)

For commercial use or if the AGPL-3.0 restrictions are not suitable for your use case, please contact us for a commercial license: zhouw3@chop.edu

## Citation

If you use MethScope, please cite (coming soon):

Fu H*, Xu H*, Lee CN, Cloud C, Deng Y, Zhou W.
**MethScope: Ultra-Fast Analysis of Sparse DNA Methylome via Recurrent Pattern Encoding**.
