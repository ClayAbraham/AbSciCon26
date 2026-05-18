param(
    [string]$InputRoot = ".\analysis\Metabolic_Profiling_v2_1\output",
    [string]$OutputDir = ".\analysis\Metabolic_Profiling_v2_1\output\cross_sample_analysis",
    [int]$MinPathwayAbundanceCount = 10,
    [int]$MinKoCount = 20,
    [int]$MinObservedKos = 3,
    [int]$MinTotalKos = 20
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

py -3 .\analysis\Metabolic_Profiling_v2_1\cross_sample_analysis.py `
    --input-root $resolvedInputRoot `
    --output-dir $resolvedOutputDir `
    --min-pathway-abundance-count $MinPathwayAbundanceCount `
    --min-ko-count $MinKoCount `
    --min-observed-kos $MinObservedKos `
    --min-total-kos $MinTotalKos
if ($LASTEXITCODE -ne 0) {
    throw "Cross-sample analysis failed"
}

Write-Host "Cross-sample analysis written to $resolvedOutputDir"
