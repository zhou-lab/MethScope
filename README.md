<h1><img src="MethScope.png" alt="MethScope logo" width="44"> MethScope</h1>

**MethScope** is an R package for ultra-fast analysis of sparse DNA methylome data using **Most Recurrent Methylation Patterns (MRMPs)**.

It supports downstream analysis for:
- Cell type annotation
- Cell type deconvolution
- Unsupervised clustering
- Cancer cell-of-origin prediction
- Missing value imputation

## Why MethScope?

Sparse single-cell and spatial methylome data are difficult to analyze directly. MethScope compresses methylation signals into MRMP-based embeddings so you can run robust and scalable downstream tasks with standard analysis workflows.

## Installation

Install from CRAN:

```r
install.packages("MethScope")
```

Or install the development version from GitHub:

```r
remotes::install_github("zhou-lab/MethScope")
```

## Quick start

```r
library(MethScope)

# 1) Generate MRMP embedding from your .cg file and MRMP reference
example_file <- "example.cg"
reference_pattern <- "Liu2021_MouseBrain.cm"
input_pattern <- GenerateInput(example_file, reference_pattern)

# 2) Predict cell types with a built-in model
pred <- PredictCellType(MethScope:::Liu2021_MouseBrain_P1000, input_pattern)

# 3) Visualize prediction results
PlotUMAP(input_pattern, pred)
```

## Tutorials and documentation

- Documentation website: [zhou-lab.github.io/MethScope](https://zhou-lab.github.io/MethScope/)
- End-to-end tutorial: [MethScope-Tutorial](https://zhou-lab.github.io/MethScope/articles/MethScope-Tutorial.html)
- Building MRMP references: [MethScope-MRMP](https://zhou-lab.github.io/MethScope/articles/MethScope-MRMP.html)

## Data resources

- Example and reference data: [zhou-lab/methscope_data](https://github.com/zhou-lab/methscope_data)
- `.cg` generation and preprocessing: [YAME](https://zhou-lab.github.io/YAME/)
- Pattern interpretation: [knowYourCG](https://www.bioconductor.org/packages/release/bioc/html/knowYourCG.html)

## Citation

If you use MethScope, please cite (comming soon):

Fu H, Xu H, Lee CN, Cloud C, Deng Y, Zhou W.  
**MethScope: Ultra-Fast Analysis of Sparse DNA Methylome via Recurrent Pattern Encoding**.
