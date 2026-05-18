param(
    [string]$InputRoot = ".\analysis\Metabolic_Profiling_v2_1\output",
    [string]$CacheDir = ".\analysis\Metabolic_Profiling_v2_1\output\kegg_reference"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$resolvedInputRoot = if ([System.IO.Path]::IsPathRooted($InputRoot)) {
    $InputRoot
} else {
    Join-Path (Get-Location).Path $InputRoot
}

$resolvedCacheDir = if ([System.IO.Path]::IsPathRooted($CacheDir)) {
    $CacheDir
} else {
    Join-Path (Get-Location).Path $CacheDir
}

[System.IO.Directory]::CreateDirectory($resolvedCacheDir) | Out-Null

py -3 .\analysis\Metabolic_Profiling_v2_1\generate_sample_pathway_ko_tables.py `
    --input-root $resolvedInputRoot `
    --cache-dir $resolvedCacheDir
if ($LASTEXITCODE -ne 0) {
    throw "Pathway KO table generation failed"
}

Write-Host "Pathway KO tables written under $resolvedInputRoot"
