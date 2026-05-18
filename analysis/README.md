# Analysis Workspaces

This directory contains downstream analysis workflows that start from compiled pathway outputs rather than raw contigs.

- `Metabolic_Profiling_v2_1/`: current Excel-first workflow for normalized matrices, charts, cross-sample comparisons, and pathway-KO summaries.
- `Metabolic_Profiling_v2_2/`: focused downstream scripts for filtered heatmaps and psychrophilic/anaerobic pathway summaries.

Generated `output/` folders and large workbook files are intentionally ignored by `.gitignore`. For reproducible releases, publish large inputs and generated outputs separately or attach them to a versioned archive.
