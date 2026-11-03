<#
.SYNOPSIS
Runs the scTCNA web app in the foreground (Ctrl+C stops it), with the same settings as the service.

.DESCRIPTION
Starts uvicorn app.main:app from the repository root with the environment made by
setup.ps1. It listens on 127.0.0.1 only; IIS forwards the site to it.

.EXAMPLE
powershell -ExecutionPolicy Bypass -File .\run.ps1
#>
[CmdletBinding()]
param(
    # default: data\app of this repository
    [string]$DataDir = '',
    [int]$Port = 8000
)

Set-StrictMode -Version 3
$ErrorActionPreference = 'Stop'

$Root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
if (-not $DataDir) { $DataDir = Join-Path $Root 'data\app' }
$DataDir = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($DataDir).TrimEnd('\', '/')
$Python = Join-Path $Root '.venv\Scripts\python.exe'
if (-not (Test-Path $Python)) { throw "$Python not found: run setup.ps1 first" }
if (-not (Test-Path (Join-Path $DataDir 'datasets.parquet'))) {
    Write-Warning "No tables in $DataDir; every /api route will answer 503"
}

$env:SCTCNA_DATA_DIR = $DataDir
# returns memory freed by Arrow tables to the OS
$env:ARROW_DEFAULT_MEMORY_POOL = 'system'
$env:PYTHONUNBUFFERED = '1'
$env:PYTHONIOENCODING = 'utf-8'

Push-Location $Root
try {
    & $Python -m uvicorn app.main:app --host 127.0.0.1 --port $Port --workers 1
} finally {
    Pop-Location
}
exit $LASTEXITCODE
