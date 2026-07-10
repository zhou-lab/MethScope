# MethScope Workflows

## Common workflow 1: generate MRMP embeddings

Goal: convert a `.cg` file plus a reference pattern definition into an `input_pattern` matrix.

Typical example:

```r
library(MethScope)

# Run from the root of a cloned zhou-lab/MethScope repository.
example_file <- "inst/extdata/example.cg"
reference_pattern <- "inst/extdata/mm10_Liu2021.cm"
input_pattern <- GenerateInput(example_file, reference_pattern)
```

Keep `inst/extdata/example.cg.idx` next to `inst/extdata/example.cg` so sample
names and sample order are available. CRAN toy files are intentionally tiny and
are not the preferred full functional test.

Use this path when the user wants:

- feature generation
- embeddings for downstream analysis
- input for prediction, training, clustering, or deconvolution

## Common workflow 2: annotate cell types with a pre-trained model

Typical example:

```r
model <- Liu2021_MouseBrain_P1000()
prediction_result <- PredictCellType(model, input_pattern)
```

Documented built-in model examples include:

- `Liu2021_MouseBrain_P1000()`
- `Zhou2025_HumanAtlas_P1000()`

Use this path when the user wants fast annotation from an existing model.

The built-in mouse brain model expects 1000 MRMP features. Use the full GitHub
`inst/extdata/mm10_Liu2021.cm` reference for realistic testing; `toy.cm` only
has a tiny subset of patterns.

## Common workflow 3: train a custom classifier

Typical example:

```r
trained_model <- Input_training(input_pattern, cell_type_label)
```

Important considerations:

- `cell_type_label` must align with rows of `input_pattern`
- The user may supply custom `xgboost` parameters
- Cross-validation may be used to search for better parameters

## Common workflow 4: visualize results

Typical examples:

```r
umap_plot <- PlotUMAP(input_pattern, prediction_result)
PlotConfusion(prediction_result, cell_type_label)
PlotF1(prediction_result, cell_type_label)
```

For the GitHub example data:

```r
example_label <- read.csv("inst/extdata/example_label.csv")
cell_type_label <- example_label$label
stopifnot(length(cell_type_label) == nrow(input_pattern))
```

The label file follows the `.cg.idx` sample order.

Use this path when the user wants QC, presentation figures, or performance summaries.

## Common workflow 5: deconvolution

Typical example:

```r
reference_input <- readRDS("inst/extdata/mm10_Liu2021_ref.rds")
input_pattern <- imputeRowMean(input_pattern)
cell_proportion <- nnls_deconv(reference_input, input_pattern)
```

Use this when the user wants bulk or mixed-sample cell type proportion estimates.
Use `hg38_Zhou2025_ref.rds` for human atlas workflows when appropriate.

## Common workflow 6: prepare `.cg` input files

MethScope uses YAME `.cg` files. The query `.cg` and reference `.cm` must be
based on the same genome build and CpG coordinate order.

Use `vignettes/MethScope-Input.Rmd` for conversion examples:

- BED-like M/U count files: align to `.cr` CpG reference, then `yame pack -f3`
- ALLC files: keep CG context, convert `mc` and `cov` to `M` and `U`, then `yame pack -f3`
- beta/fraction tables: use `yame pack -f4`
- binary tracks: use `yame pack -fb`

Point users to KYCGKB `.cr` references for genome-specific CpG order:

- hg38: `https://github.com/zhou-lab/KYCGKB_hg38`
- mm10: `https://github.com/zhou-lab/KYCGKB_mm10`

## Common workflow 7: clustering downstream of MethScope

Once `input_pattern` is available, users can cluster with external tools such as Seurat. The MethScope tutorial demonstrates this as a downstream workflow rather than a package-native clustering API.

## Common workflow 8: methscope-cli C implementation

Use the R package for interactive analysis, plotting, training, and R object
workflows. Use `methscope-cli` for standalone shell, batch, or HPC workflows
without an R runtime.

Main command patterns:

```bash
methscope predict query.cg model.ubjx > prediction.tsv
methscope matrix query.cg model_or_mrmp > input_pattern.tsv
methscope deconv mixture.cg reference.refx > deconv.tsv
methscope upscale -o reconstructed.cg upscale_model.updecx sparse_input.cg
```

Bundle formats:

- `.ubjx`: classifier bundle for `methscope predict`
- `.refx`: deconvolution reference bundle for `methscope deconv`
- `.updecx`: upscaling decoder bundle for `methscope upscale`

Self-contained bundles carry their MRMP feature definition, so users do not
need to hand-match a separate `.cm` file for those CLI workflows.

## Troubleshooting heuristics

- Missing or incompatible inputs: verify `.cg` and reference `.cm` paths first
- Poor performance on large input: recommend splitting `.cg` files and parallelizing generation
- Unexpected prediction labels: verify the model and reference are from the right organism or atlas
- Shape mismatches: check whether rows and columns are oriented as the downstream function expects
- Missing sample names or label order mismatch: verify the `.cg.idx` file is next to the `.cg` file
- Deconvolution with missing values: run `imputeRowMean(input_pattern)` first
