# MethScope

**MethScope** identifies Most Recurrent Methylation Patterns (MRMPs) as the basis to encode latent representations for interpretable analysis of DNA methylation data, in particular single cell and spatial methylome. The MRMPs embeddings will support automatic **cell annotation**, **bulk deconvolution**, **unsupervised clustering**, and **cancer cell of origin prediction**.

## System Requirements

### Software Dependencies

- **R** >= 4.0 (tested on R 4.5.2 and R 4.6.0)
- **Operating systems tested**: macOS, Linux (Ubuntu)
- R package dependencies (installed automatically): `xgboost`, `dplyr`, `tidyr`, `stringr`, `caret`, `doParallel`, `parallel`, `ggplot2`, `uwot`, `magrittr`, `FNN`, `data.table`, `nnls`

### Hardware Requirements

No non-standard hardware is required. MethScope runs on a standard laptop or desktop CPU. No GPU is needed.

## Installation

Install from CRAN:

``` r
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

For full documentation and advanced usage (deconvolution, clustering, training custom models), see the [MethScope website](https://zhou-lab.github.io/MethScope/).

## License

This project is licensed under the [GNU Affero General Public License v3.0 (AGPL-3.0)](LICENSE.md).

Copyright (c) 2025 Hongxiang Fu and Wanding Zhou (zhouw3@chop.edu)

For commercial use or if the AGPL-3.0 restrictions are not suitable for your use case, please contact us for a commercial license: zhouw3@chop.edu

## Citation

If you use MethScope, kindly cite (coming soon):

Hongxiang Fu, Chin Nien Lee, Cameron Cloud, Hao Xu, Yanxiang Deng, Wanding Zhou, MethScope: Ultra-fast Analysis of Sparse DNA Methylome via Recurrent Pattern Encoding.