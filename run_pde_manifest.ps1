# Script 100 end to end: first-pair full-grid forecasts, cell selection, later-pair selected-cell forecasts, analysis.
# 4 processes x 2 torch threads (multiprocessing.spawn is refused on this machine). Resumable: finished pairs are cached.
# Usage: powershell -ExecutionPolicy Bypass -File run_pde_manifest.ps1
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
$out = Join-Path $root "output\pde_manifest"
New-Item -ItemType Directory -Force $out | Out-Null
$n = 4

function Run-Stage($stage) {
    $procs = @()
    for ($i = 0; $i -lt $n; $i++) {
        $log = Join-Path $out "$stage-shard$i.log"
        $procs += Start-Process -FilePath "python" -ArgumentList @("-u", "src\100_pde_manifest.py", "--stage", $stage, "--shard", $i, "--n-shards", $n) `
            -RedirectStandardOutput $log -RedirectStandardError "$log.err" -WindowStyle Hidden -PassThru
    }
    $procs | Wait-Process
}

Run-Stage "forecast"
& python -u src\100_pde_manifest.py --stage select *> (Join-Path $out "select.log")
Run-Stage "selected"
& python -u src\100_pde_manifest.py --stage analyze *> (Join-Path $out "analyze.log")
"done" | Out-File (Join-Path $out "DONE.txt")
