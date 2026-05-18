param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$PipelineArgs
)

$ErrorActionPreference = "Stop"

$ProjectRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
Set-Location $ProjectRoot

$EnvName = (Get-Content (Join-Path $ProjectRoot ".conda-env") -Raw).Trim()
if (-not $EnvName) {
    throw "Missing conda environment name in .conda-env"
}

& conda run -n $EnvName python (Join-Path $ProjectRoot "main.py") @PipelineArgs
exit $LASTEXITCODE
