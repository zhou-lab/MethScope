# MethScope Core Instructions

Use this skill when the user is working with the MethScope R package, its tutorial, vignettes, reference data, or downstream analysis workflows.

## Package identity

- Package name: `MethScope`
- Repository root: the directory containing `DESCRIPTION`
- Tutorial site: `https://zhou-lab.github.io/MethScope/`
- Main vignette: `vignettes/MethScope-Tutorial.Rmd`

## What MethScope is for

MethScope identifies Most Recurrent Methylation Patterns (MRMPs) and uses those embeddings for sparse DNA methylome analysis, especially single-cell and spatial methylome data.

Primary use cases:

- Generate cell, pixel, or sample by MRMP embedding matrices
- Perform cell type annotation with pre-trained models
- Train custom classification models
- Visualize prediction and embedding results
- Estimate cell type proportions by deconvolution
- Support clustering and exploratory downstream analysis

## Trigger conditions

Use this skill if the request involves any of the following:

- MethScope package functions or internals
- `.cg`, `.cg.idx`, or `.cm` files used by MethScope
- MRMP embeddings
- Cell type prediction with `PredictCellType`
- MethScope model training with `Input_training`
- Deconvolution with `nnls_deconv`
- MethScope vignettes, examples, or pkgdown documentation

## Default working assumptions

- Prefer the package's own functions and vignettes over inventing a custom workflow
- Prefer examples from the package tutorial when the user asks for usage help
- If a file format detail is unclear, point back to the documented source rather than guessing
- Keep the distinction clear between reference pattern files (`.cm`) and methylation input files (`.cg`)

## Canonical first questions to answer internally

Before proposing commands or code, identify:

1. Is the user trying to generate embeddings, annotate cells, train a model, deconvolve, or visualize?
2. Do they already have a `.cg` input file and a reference `.cm` pattern file?
3. Are they using a built-in pre-trained model or training a custom model?
4. Do they want package usage help, code changes inside the package, or interpretation of outputs?

## Preferred information sources

Consult these in roughly this order:

1. The relevant vignette or pkgdown article
2. The exported function docs in `man/` or `docs/reference/`
3. The implementation in `R/`
4. Example analysis scripts in `analysis/`

## Guardrails

- Do not assume all methylation formats are interchangeable with MethScope inputs
- Do not hallucinate built-in models or reference datasets beyond those documented
- Do not recommend changing package internals when the request is only about package usage
- If an answer depends on local files or model artifacts not present, say so clearly
