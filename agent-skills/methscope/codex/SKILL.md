---
name: methscope
description: Use when working with the MethScope R package, MethScope documentation, MRMP embedding workflows, methylation input files, cell type prediction, training, deconvolution, visualization, the methscope-cli C implementation, pretrained bundles, or package internals.
---

# MethScope Skill

Use this skill for MethScope package usage, analysis, documentation, and development tasks.

Read these files as needed:

- `../core/INSTRUCTIONS.md` for scope, triggers, and guardrails
- `../core/WORKFLOWS.md` for common package workflows
- `../core/REFERENCES.md` for function names, files, and reference paths

## Operating guidance

- Prefer documented MethScope workflows over generic methylation-analysis advice
- Use the tutorial and vignette examples as the default source of usage patterns
- Keep file format assumptions explicit, especially for `.cg` and `.cm`
- Keep the GitHub full example data distinct from the tiny CRAN toy files
- Treat the repository root as the directory containing `DESCRIPTION`
- For code changes, inspect the relevant files under `R/`, `vignettes/`, `man/`, `docs/`, and `inst/extdata/`
