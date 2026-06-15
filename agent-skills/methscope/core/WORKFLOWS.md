# MethScope Workflows

## Common workflow 1: generate MRMP embeddings

Goal: convert a `.cg` file plus a reference pattern definition into an `input_pattern` matrix.

Typical example:

```r
library(MethScope)

example_file <- "example.cg"
reference_pattern <- "mm10_Liu2021.cm"
input_pattern <- GenerateInput(example_file, reference_pattern)
```

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

Use this path when the user wants QC, presentation figures, or performance summaries.

## Common workflow 5: deconvolution

Typical example:

```r
reference_input <- readRDS("2021Liu_reference_pattern.rds")
cell_proportion <- nnls_deconv(reference_input, input_pattern)
```

Use this when the user wants bulk or mixed-sample cell type proportion estimates.

## Common workflow 6: clustering downstream of MethScope

Once `input_pattern` is available, users can cluster with external tools such as Seurat. The MethScope tutorial demonstrates this as a downstream workflow rather than a package-native clustering API.

## Troubleshooting heuristics

- Missing or incompatible inputs: verify `.cg` and reference `.cm` paths first
- Poor performance on large input: recommend splitting `.cg` files and parallelizing generation
- Unexpected prediction labels: verify the model and reference are from the right organism or atlas
- Shape mismatches: check whether rows and columns are oriented as the downstream function expects
