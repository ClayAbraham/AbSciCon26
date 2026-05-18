# Metagenomic Metabolic Profiling

Python workflows for metagenomic metabolic profiling, from contig-level annotation through pathway and KO-level downstream analysis.

This repository contains two active workflow layers:

- A sequence-first V1 pipeline that starts from contig FASTA files, predicts proteins, assigns KEGG Orthology identifiers, and scores pathway completeness.
- Excel-first downstream analysis workspaces for compiled pathway outputs, cross-sample comparisons, focused heatmaps, and pathway/gene summaries.

## Repository Layout

| Path | Purpose |
| --- | --- |
| `metabolic_pipeline/` | Primary V1 Python package for contig-to-pathway profiling. |
| `analysis/Metabolic_Profiling_v2_1/` | Current Excel-first analysis workflow. |
| `analysis/Metabolic_Profiling_v2_2/` | Focused downstream scripts for filtered pathway, psychrophilic, and anaerobic analyses. |
| `config/` | Reusable configuration templates. |
| `scripts/` | Local and cluster launch scripts. |
| `docs/` | Lightweight static documentation / GitHub Pages files. |
| `tests/` | Basic smoke tests. |
| `archive/` | Earlier workflow material, command logs, and provenance notes. |

See `REPOSITORY_GUIDE.md` for repository hygiene and release guidance.

## V1 Sequence Pipeline

For each input sample, the V1 pipeline:

1. Reads contig FASTA files (`.fasta`, `.fa`, or `.fna`).
2. Runs Prodigal to predict proteins.
3. Runs KOfamScan on predicted proteins, or reads precomputed KOfam output.
4. Optionally parses precomputed eggNOG annotations.
5. Merges annotations and summarizes KO counts.
6. Scores user-defined pathways.
7. Writes tabular reports and run metadata.

### Environment

The expected conda environment name is stored in `.conda-env`:

```text
metabolic-pipeline
```

Create the environment with:

```bash
conda env create -f environment.yml
```

Required external tools:

- Prodigal
- KOfamScan and its database files, if `run_kofam: true`

`environment.yml` installs Prodigal. KOfamScan database paths should be supplied in the config file.

### Local Run

Linux/macOS:

```bash
./scripts/run_pipeline.sh \
  --input /path/to/contigs.fasta \
  --output results \
  --pathways /path/to/pathways.json \
  --run-kofam \
  --overwrite
```

Windows PowerShell:

```powershell
.\scripts\run_pipeline.ps1 `
  --input C:\path\to\contigs.fasta `
  --output results `
  --pathways C:\path\to\pathways.json `
  --run-kofam `
  --overwrite
```

If `--run-kofam` is omitted, expected precomputed KOfam files must already exist at:

```text
results/kofam/<sample_id>.kofam.tsv
```

### Config-Driven Run

Copy and edit the example config:

```bash
cp config/config.cluster.example.yaml config/config.cluster.yaml
```

Run with:

```bash
./scripts/run_pipeline.sh --config config/config.cluster.yaml
```

Important config fields:

- `input_path`: contig FASTA file or directory
- `output_dir`: output directory
- `required_conda_env`: expected conda environment name
- `prodigal_executable`: Prodigal executable name or full path
- `prodigal_mode`: `meta` or `single`
- `pathway_definition_file`: pathway JSON
- `kofamscan_executable`: KOfamScan executable path
- `kofam_profiles_dir`: KOfam profiles directory
- `kofam_ko_list`: KOfam `ko_list` path

### Cluster Run

Edit `scripts/submit_metabolic_pipeline.sbatch` for the cluster project path, then submit:

```bash
sbatch scripts/submit_metabolic_pipeline.sbatch
```

## V2.1 Excel Analysis

Use `analysis/Metabolic_Profiling_v2_1/` when starting from compiled Excel pathway workbooks rather than raw sequence data.

```powershell
conda env create -f .\analysis\Metabolic_Profiling_v2_1\environment.v2.yml
conda activate metabolic-profiling-v2
py -3 .\analysis\Metabolic_Profiling_v2_1\run_analysis.py --config .\analysis\Metabolic_Profiling_v2_1\config.example.yaml
```

This workflow can normalize long- or wide-format worksheets, compute feature and sample statistics, compare samples, render charts, and run threshold-aware cross-sample analyses.

## Outputs

Typical V1 outputs include:

- `prodigal/<sample_id>.faa`
- `prodigal/<sample_id>.gff`
- `kofam/<sample_id>.kofam.tsv`
- `annotations_merged.tsv`
- `ko_summary.tsv`
- `pathway_summary.tsv`
- `pathway_long.tsv`
- `run.log`
- `config_snapshot.yaml`
- `versions.txt`

Generated output folders, local environments, bytecode caches, and large workbook/chart files are ignored by `.gitignore` by default. Publish large generated artifacts separately when creating a reproducible release.

## Testing

Basic syntax/import checks:

```bash
python -m compileall metabolic_pipeline analysis/Metabolic_Profiling_v2_1 analysis/Metabolic_Profiling_v2_2 tests
```

The V1 CLI requires the project dependencies from `environment.yml`.

## Related Resources

Add project-specific links here before publication:

- People and affiliations: `PEOPLE_AND_AFFILIATIONS.md`
- Project website:
- Dataset archive:
- Paper or preprint:
- Documentation:

## Citation

If you use this project in academic work, cite the repository metadata in `CITATION.cff`. Replace the placeholder `repository-code` URL with the final GitHub repository URL before release.

## License

This project is distributed under the MIT license. Confirm that this is the intended license before public publication, especially for institutional or collaborative work.
