# Metabolic Profiling V2.1

This folder is a self-contained analysis workspace for Excel-based metabolic profiling outputs.

V2.1 is the threshold-aware branch of V2. It is intended for stricter cross-sample analysis where low-signal mapping values can be filtered before clustering and sample-to-sample comparison.

The goal of V2.1 is different from the original pipeline:

- V1 starts from sequence data and annotation tools.
- V2.1 starts from compiled Excel workbooks and produces charts, graphs, and statistical summaries.

## What this version can do

Given one or more worksheet tables in an Excel workbook, V2.1 can:

- normalize long-format or wide-format pathway tables into sample-by-feature matrices
- calculate feature-level statistics
- calculate sample-level statistics
- calculate pathway prevalence across samples
- calculate sample-to-sample correlation and distance matrices
- run PCA on samples when there is enough variance
- run optional group-based tests if metadata is supplied
- export charts for abundance and completeness-style datasets
- apply lower-value thresholds during cross-sample analysis to suppress likely mapping noise

## Recommended workbook structure

This script is designed to support either of these layouts:

### Long format

| pathway_name | sample_id | abundance |
| --- | --- | --- |
| Glycolysis / Gluconeogenesis | 1A | 6664 |
| Glycolysis / Gluconeogenesis | 5A | 7122 |

### Wide format

| pathway_name | 1A | 5A | 14A |
| --- | --- | --- | --- |
| Glycolysis / Gluconeogenesis | 6664 | 7122 | 5901 |

For pathway completeness, the value column might be something like `completeness_pct`.

## Folder layout

- `run_analysis.py` entry point
- `config.example.yaml` analysis template
- `input/` for your Excel workbook
- `command_logs/` for your saved workflow notes
- `output/` for generated tables and figures

## Install analysis dependencies

Create or update an environment with the packages listed in `requirements.txt`, or create a dedicated conda environment from `environment.v2.yml`.

```powershell
conda env create -f .\analysis\Metabolic_Profiling_v2_1\environment.v2.yml
conda activate metabolic-profiling-v2
```

Then run:

```powershell
py -3 .\analysis\Metabolic_Profiling_v2_1\run_analysis.py --config .\analysis\Metabolic_Profiling_v2_1\config.example.yaml
```

For the current sample workbook already placed in this folder, use:

```powershell
py -3 .\analysis\Metabolic_Profiling_v2_1\run_analysis.py --config .\analysis\Metabolic_Profiling_v2_1\config.1A.yaml
```

## Batch runs

Render all per-sample charts:

```powershell
powershell -ExecutionPolicy Bypass -File .\analysis\Metabolic_Profiling_v2_1\generate_all_charts.ps1
```

Run the threshold-aware cross-sample analysis:

```powershell
powershell -ExecutionPolicy Bypass -File .\analysis\Metabolic_Profiling_v2_1\generate_cross_sample_analysis.ps1
```

Current default cross-sample thresholds in V2.1:

- pathway abundance values below `10` are set to zero
- KO count values below `20` are set to zero
- pathway completeness entries are removed when `total_kos < 20`
- pathway completeness entries are removed when estimated `observed_kos < 3`

## Next step

For a new analysis, place these files here:

- your Excel workbook in `analysis/Metabolic_Profiling_v2_1/input/`
- your command log copy in `analysis/Metabolic_Profiling_v2_1/command_logs/`

Then we can tune the config to the real sheet names and columns.

Note:

- the current `1A_compiled_pathway_data.xlsx` workbook is already configured through `config.1A.yaml`
- the current command log has also been recognized in the V2.1 folder
