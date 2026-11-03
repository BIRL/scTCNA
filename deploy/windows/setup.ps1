<#
.SYNOPSIS
Installs the scTCNA web app on Windows.

.DESCRIPTION
Uses .venv in the repository (creates it with Python 3.12 if it does not exist),
installs requirements-app.txt and checks that the tables in -DataDir load. The app
runs from the repository checkout; IIS is not touched. Safe to run again for
updates; stop the scTCNA service first.

.EXAMPLE
powershell -ExecutionPolicy Bypass -File .\setup.ps1
#>
[CmdletBinding()]
param(
    # default: data\app of this repository
    [string]$DataDir = '',
    # python.exe for a new environment; default: Python 3.12 from the py launcher
    [string]$Python = ''
)

Set-StrictMode -Version 3
$ErrorActionPreference = 'Stop'

# Windows PowerShell does not stop on a failing native command, so check exit codes.
function Invoke-Checked {
    param([Parameter(Mandatory = $true)][string]$FilePath, [string[]]$ArgumentList = @())
    & $FilePath @ArgumentList
    if ($LASTEXITCODE -ne 0) { throw "$FilePath $($ArgumentList -join ' ') failed with exit code $LASTEXITCODE" }
}

# The pinned numpy needs Python 3.12 or newer; the app is built and tested with 64-bit Python 3.12 (Dockerfile).
function Assert-SupportedPython {
    param([Parameter(Mandatory = $true)][string]$Exe, [Parameter(Mandatory = $true)][string]$Hint)
    $v = (& $Exe -c 'import sys; print(sys.version_info[0] * 100 + sys.version_info[1], 64 if sys.maxsize > 2**32 else 32)').Split(' ')
    if ($v[0] -ne '312' -or $v[1] -ne '64') { throw "$Exe is Python $($v[0]) ($($v[1])-bit); scTCNA needs 64-bit Python 3.12. $Hint" }
}

$Svc = Get-Service scTCNA -ErrorAction SilentlyContinue
if ($Svc -and $Svc.Status -ne 'Stopped') { throw 'Stop the service first: Stop-Service scTCNA' }

$Root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$Venv = Join-Path $Root '.venv'
$VenvPython = Join-Path $Venv 'Scripts\python.exe'
if (-not $DataDir) { $DataDir = Join-Path $Root 'data\app' }
$DataDir = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($DataDir).TrimEnd('\', '/')

Write-Host '== Python environment'
if (-not (Test-Path $VenvPython)) {
    if (-not $Python) {
        if (-not (Get-Command py -ErrorAction SilentlyContinue)) { throw 'The py launcher was not found: install Python 3.12 (64-bit) from python.org, or pass -Python <path to python.exe>' }
        $Python = (& py -3.12 -c 'import sys; print(sys.executable)')
        if ($LASTEXITCODE -ne 0) { throw 'Python 3.12 is not installed: install it from python.org (64-bit), or pass -Python <path to python.exe>' }
    }
    Assert-SupportedPython $Python 'Pass -Python <path to python.exe> of a supported version.'
    Invoke-Checked $Python @('-m', 'venv', $Venv)
}
Assert-SupportedPython $VenvPython 'Rename .venv to .venv.old and rerun setup.ps1 to create a new one.'
Invoke-Checked $VenvPython @('-m', 'pip', 'install', '--disable-pip-version-check', '-r', (Join-Path $Root 'requirements-app.txt'))

Write-Host "== Tables in $DataDir"
# sctcna is imported from the checkout, as when the app runs, so run from the repository root.
Push-Location $Root
try {
    Invoke-Checked $VenvPython @('-c', "import sys; from pathlib import Path; from sctcna.web import load_store; s, e = load_store(Path(sys.argv[1])); sys.exit(e) if s is None else print('tables load:', s.health()['n_datasets'], 'datasets')", $DataDir)
} finally {
    Pop-Location
}

Write-Host ''
Write-Host 'Done. Test the app in the foreground with:'
Write-Host "  powershell -ExecutionPolicy Bypass -File `"$PSScriptRoot\run.ps1`" -DataDir `"$DataDir`""
Write-Host 'then install the service and the IIS site as described in README_WINDOWS.md.'
