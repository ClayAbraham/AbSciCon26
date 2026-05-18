# Repository Guide

This repository contains three related workflow layers.

## Recommended public repository contents

Track source code, small example inputs, documentation, environment files, and reproducibility scripts. Do not track local virtual environments, Python bytecode, large Excel workbooks, generated charts, or routine output folders unless a specific release requires those artifacts.

## Workflow layers

- `metabolic_pipeline/`: V1 sequence-first pipeline. Starts from contig FASTA files, predicts proteins with Prodigal, assigns KOs with KOfamScan or precomputed KOfam output, and scores pathway completeness.
- `analysis/Metabolic_Profiling_v2_1/`: current Excel-first analysis workspace. Starts from compiled pathway workbooks and produces normalized matrices, charts, cross-sample statistics, and pathway-KO summaries.
- `analysis/Metabolic_Profiling_v2_2/`: downstream focused analyses for filtered pathway heatmaps and psychrophilic/anaerobic pathway summaries.

`archive/Metabolic_Profiling_v2/` is retained for provenance as an earlier V2 workspace. Prefer `analysis/Metabolic_Profiling_v2_1/` for new Excel-based runs.

## Before first GitHub upload

1. Replace the placeholder `repository-code` URL in `CITATION.cff`.
2. Confirm the MIT license is the intended license for the project.
3. Decide whether generated `output/` artifacts should be excluded from the first commit or attached separately as a release/archive.
4. Check that no private data, unpublished sample metadata, credentials, or local machine paths are included.
