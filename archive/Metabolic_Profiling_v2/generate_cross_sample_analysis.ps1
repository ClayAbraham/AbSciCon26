param(
    [string]$InputRoot = ".\Metabolic_Profiling_v2\output",
    [string]$OutputDir = ".\Metabolic_Profiling_v2\output\cross_sample_analysis"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$resolvedInputRoot = if ([System.IO.Path]::IsPathRooted($InputRoot)) {
    $InputRoot
} else {
    Join-Path (Get-Location).Path $InputRoot
}

$resolvedOutputDir = if ([System.IO.Path]::IsPathRooted($OutputDir)) {
    $OutputDir
} else {
    Join-Path (Get-Location).Path $OutputDir
}

[System.IO.Directory]::CreateDirectory($resolvedOutputDir) | Out-Null

py -3 .\Metabolic_Profiling_v2\cross_sample_analysis.py --input-root $resolvedInputRoot --output-dir $resolvedOutputDir
if ($LASTEXITCODE -ne 0) {
    throw "Cross-sample analysis failed"
}

Write-Host "Cross-sample analysis written to $resolvedOutputDir"
