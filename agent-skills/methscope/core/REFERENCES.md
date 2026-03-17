# MethScope References

## Important local files

- Package root: repository root containing `DESCRIPTION`
- Description: `DESCRIPTION`
- README: `README.md`
- Tutorial vignette: `vignettes/MethScope-Tutorial.Rmd`
- Pkgdown tutorial page: `docs/articles/MethScope-Tutorial.html`

## Exported functions

- `GenerateInput`
- `GenerateInput_withCoverage`
- `GenerateReference`
- `Input_training`
- `PredictCellType`
- `PlotUMAP`
- `PlotUMAP_fixedwindow`
- `PlotConfusion`
- `PlotF1`
- `nnls_deconv`
- `confidence_score`
- `confidence_score_top95`
- `filter_cell`
- `imputeRowMean`
- `smooth_matrix`

## Relevant source files

- `R/GenerateInput.R`
- `R/PredictCellType.R`
- `R/ModelTraining.R`
- `R/VisualizeOutput.R`
- `R/utils.R`

## Built-in example and reference data

- `inst/extdata/example.cg`
- `inst/extdata/example.cg.idx`
- `inst/extdata/toy.cg`
- `inst/extdata/toy.cg.idx`
- `inst/extdata/toy.cm`
- `inst/extdata/Liu2021_MouseBrain.cm`
- `inst/extdata/Loyfer2023_HumanAtlas.cm`
- `inst/extdata/Zhou2025_HumanAtlas.cm`

## External references

- Tutorial site: `https://zhou-lab.github.io/MethScope/`
- Data repository: `https://github.com/zhou-lab/methscope_data`
- YAME documentation for `.cg` generation and splitting: `https://zhou-lab.github.io/YAME/`
- knowYourCG package page: `https://www.bioconductor.org/packages/release/bioc/html/knowYourCG.html`

## Agent behavior notes

- When giving usage help, start from the tutorial examples
- When answering implementation questions, inspect `R/` files directly
- When asked to reproduce or debug a workflow, prefer package example data before using user data
