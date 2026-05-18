param(
    [string]$InputDir = ".\Metabolic_Profiling_v2\input",
    [string]$OutputBaseDir = ".\Metabolic_Profiling_v2\output"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$resolvedInputDir = (Resolve-Path $InputDir).Path
$resolvedOutputBaseDir = [System.IO.Path]::GetFullPath((Join-Path (Get-Location) $OutputBaseDir))
[System.IO.Directory]::CreateDirectory($resolvedOutputBaseDir) | Out-Null

$workbooks = Get-ChildItem -Path $resolvedInputDir -Filter *.xlsx | Sort-Object Name
if (-not $workbooks) {
    throw "No .xlsx workbooks were found in $resolvedInputDir"
}

foreach ($workbook in $workbooks) {
    $sampleId = $workbook.BaseName.Split('_')[0]
    $sampleOutputDir = Join-Path $resolvedOutputBaseDir ($sampleId + "_charts")
    & powershell -ExecutionPolicy Bypass -File .\Metabolic_Profiling_v2\generate_1A_charts.ps1 `
        -WorkbookPath $workbook.FullName `
        -OutputDir $sampleOutputDir `
        -SampleId $sampleId
    if ($LASTEXITCODE -ne 0) {
        throw "Chart generation failed for sample $sampleId"
    }
}

Write-Host "Processed $($workbooks.Count) workbooks into $resolvedOutputBaseDir"
