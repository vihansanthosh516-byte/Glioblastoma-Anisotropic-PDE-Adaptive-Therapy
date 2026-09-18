#!/usr/bin/env pwsh
# run_all.ps1 — PowerShell equivalent of run_all.sh (Month 7 -> Month 10).
# Chains the PDE cohort scripts so that Month 10's ingest step always
# finds fresh metric JSONs.
#
# Usage:
#   powershell.exe -File run_all.ps1           # run the full chain
#   powershell.exe -File run_all.ps1 -Month10  # only run script 45 (assumes 42-44 already done)
# Notes:
#   - Requires a Python venv at ./venv or via $PYTHON env var.
#   - Set GBM_PROJECT_ROOT if not using default path.

# Resolve project root (directory containing this script).
$scriptDir = $MyInvocation.MyCommand.Definition | Split-Path -Parent
$scriptDir = Get-Item -Path $scriptDir -ErrorAction SilentlyContinue
if (-not $scriptDir) {
    $scriptDir = Get-Location
}
$projectRoot = $scriptDir.Parent.Parent
if (-not $projectRoot) {
    Write-Error "Could not determine project root." ; exit 1
}
cd $projectRoot

# Determine Python interpreter.
$py = $env:PYTHON
if (-not $py -or (-not (Test-Path "$py"))) {
    $py = "$env:HOMEPATH\20206 science fair\venv\Scripts\python.exe"
    if (-not (Test-Path $py)) {
        Write-Error "Python interpreter not found. Set PYTHON env var or activate a venv." ; exit 1
    }
}

# Define scripts and months.
$scripts = @(
    "42_anisotropic_pde",
    "43_stromal_feedback",
    "44_adaptive_therapy",
    "46_sensitivity_analysis",
    "47_optimal_control",
    "48_3d_extension",
    "45_validation_synthesis"
)

$months = @{
    "42_anisotropic_pde"      = "7"
    "43_stromal_feedback"     = "8"
    "44_adaptive_therapy"     = "9"
    "46_sensitivity_analysis" = "Phase 2b"
    "47_optimal_control"      = "Phase 3"
    "48_3d_extension"         = "Phase 3D"
    "45_validation_synthesis" = "10"
}

# --month10 switch => only run script 45
if ($args -contains "-Month10") -or ($args -contains "-month10") -or ($args -contains "/Month10") -or ($args -contains "/month10") {
    $scripts = @("45_validation_synthesis")
}

Write-Host "==================================================================" -ForegroundColor Cyan
Write-Host " RUNNING COHORT PIPELINE  (Month 7 -> Month 10)" -ForegroundColor Cyan
Write-Host "==================================================================" -ForegroundColor Cyan
Write-Host "  Interpreter : $py"
Write-Host "  Project root: $projectRoot"
Write-Host "  Scripts     : $($scripts -join ', ')"
Write-Host "==================================================================" -ForegroundColor Cyan

foreach ($s in $scripts) {
    src = "src/${s}.py"
    if (-not (Test-Path $src)) {
        Write-Host "[run_all] ERROR: missing $src" -ForegroundColor Red ; exit 1
    }
    Write-Host "########## MONTH $($months[$s]) — $s ##########"
    & $py "$src"
    if ($LASTEXITCODE -ne 0) {
        Write-Host "[run_all] FAILED at $src" -ForegroundColor Red ; exit 1
    }
}

Write-Host "`n"
Write-Host "==================================================================" -ForegroundColor Cyan
Write-Host " ALL SCRIPTS COMPLETED SUCCESSFULLY" -ForegroundColor Cyan
Write-Host "==================================================================" -ForegroundColor Cyan
Write-Host " Month 10 deliverables:"
Write-Host "   output/master_cohort_summary.json"
Write-Host "   output/master_cohort_synthesis.png"
Write-Host "   output/POSTER_KEY_FINDINGS.md"
Write-Host "   output/MONTH10_AUDIT.md"
Write-Host "   output/isotropic_baseline_metrics.json"
Write-Host "==================================================================" -ForegroundColor Cyan