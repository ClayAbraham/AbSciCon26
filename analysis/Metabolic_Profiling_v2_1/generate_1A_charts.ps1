param(
    [string]$WorkbookPath = ".\analysis\Metabolic_Profiling_v2_1\input\1A_compiled_pathway_data.xlsx",
    [string]$OutputDir = ".\analysis\Metabolic_Profiling_v2_1\output\1A_charts",
    [string]$SampleId = "1A"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$resolvedWorkbookPath = (Resolve-Path $WorkbookPath).Path
$resolvedOutputDir = if ([System.IO.Path]::IsPathRooted($OutputDir)) {
    $OutputDir
} else {
    Join-Path (Get-Location).Path $OutputDir
}
[System.IO.Directory]::CreateDirectory($resolvedOutputDir) | Out-Null

py -3 .\analysis\Metabolic_Profiling_v2_1\extract_1A_data.py --workbook $resolvedWorkbookPath --output-dir $resolvedOutputDir
if ($LASTEXITCODE -ne 0) {
    throw "Extraction failed for sample $SampleId"
}
py -3 .\analysis\Metabolic_Profiling_v2_1\render_1A_charts.py --input-dir $resolvedOutputDir --output-dir $resolvedOutputDir --sample-id $SampleId
if ($LASTEXITCODE -ne 0) {
    throw "Rendering failed for sample $SampleId"
}

Write-Host "Charts and summary files written to $resolvedOutputDir"
