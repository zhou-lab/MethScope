# MethScope References

## Important local files

- Package root: repository root containing `DESCRIPTION`
- Description: `DESCRIPTION`
- README: `README.md`
- Tutorial vignette: `vignettes/MethScope-Tutorial.Rmd`
- Input preparation vignette: `vignettes/MethScope-Input.Rmd`
- MRMP reference vignette: `vignettes/MethScope-MRMP.Rmd`
- C CLI vignette: `vignettes/methscope-cli.Rmd`
- Pretrained models vignette: `vignettes/pretrained-models.Rmd`
- Pkgdown tutorial page: `docs/articles/MethScope-Tutorial.html`
- Pkgdown input page: `docs/articles/MethScope-Input.html`
- Pkgdown CLI page: `docs/articles/methscope-cli.html`
- Pkgdown pretrained models page: `docs/articles/pretrained-models.html`

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
- `inst/extdata/mm10_Liu2021.cm`
- `inst/extdata/hg38_Loyfer2023.cm`
- `inst/extdata/hg38_Zhou2025.cm`
- `inst/extdata/example_label.csv`
- `inst/extdata/mm10_Liu2021_ref.rds`
- `inst/extdata/hg38_Zhou2025_ref.rds`

Notes:

- `inst/extdata/example.cg` is the full GitHub example for functional testing.
- CRAN size limits mean the CRAN package contains tiny `toy.cg` and `toy.cm` files for minimal validation, not realistic prediction tests.
- `example.cg.idx` should be kept in the same folder as `example.cg`; `example_label.csv` follows the `.cg.idx` sample order.
- GitHub `.cm` files are named by genome build and source dataset: `mm10_Liu2021.cm`, `hg38_Zhou2025.cm`, and `hg38_Loyfer2023.cm`.
- Built-in mouse brain prediction uses the first 1000 patterns from the full `mm10_Liu2021.cm` reference.

## External references

- Tutorial site: `https://zhou-lab.github.io/MethScope/`
- MethScope input tutorial: `https://zhou-lab.github.io/MethScope/articles/MethScope-Input.html`
- methscope-cli tutorial: `https://zhou-lab.github.io/MethScope/articles/methscope-cli.html`
- Pretrained models and bundles: `https://zhou-lab.github.io/MethScope/articles/pretrained-models.html`
- methscope-cli repository: `https://github.com/zhou-lab/methscope-cli`
- Data and model bundle repository: `https://github.com/zhou-lab/methscope_data`
- YAME documentation for `.cg` generation and splitting: `https://zhou-lab.github.io/YAME/`
- YAME storage guide: `https://zhou-lab.github.io/YAME/docs/storage.html`
- ALLC input format reference: `https://lhqing.github.io/ALLCools/start/input_files.html`
- KYCGKB hg38 `.cr` references: `https://github.com/zhou-lab/KYCGKB_hg38`
- KYCGKB mm10 `.cr` references: `https://github.com/zhou-lab/KYCGKB_mm10`
- knowYourCG package page: `https://www.bioconductor.org/packages/release/bioc/html/knowYourCG.html`

## Agent behavior notes

- When giving usage help, start from the tutorial examples
- When answering implementation questions, inspect `R/` files directly
- When asked to reproduce or debug a workflow, prefer GitHub `example.cg` before user data, unless the task is explicitly about CRAN installation sanity checks
- For deconvolution examples, prefer bundled `mm10_Liu2021_ref.rds` or `hg38_Zhou2025_ref.rds` instead of external loose files
- For CLI questions, distinguish R package functions from `methscope-cli` commands and bundle formats
