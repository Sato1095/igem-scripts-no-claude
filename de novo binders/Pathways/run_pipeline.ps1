# run_pipeline.ps1 -- regenerate every artifact of the Boolean model project.
#
#   .\run_pipeline.ps1
#   .\run_pipeline.ps1 -Refetch          # re-download from OmniPath (not advised)
#
# Requires: Python 3.8+ (standard library only; no packages to install).
#
# Inputs that must already be present in this folder:
#   model_scope_v1.txt        declared node scope
#   decisions_final.csv       25 curation decisions
#   rules_v3_unfitted.txt     Boolean rules
#   screen_criteria_v2.md     pre-registered criteria (read by you, not by code)
#   run01\raw\                archived OmniPath snapshot  <-- do not delete
#
# Everything else in this folder is regenerated.

param(
    [string]$Python = "C:\Users\dimko\AppData\Local\Programs\Python\Python314\python.exe",
    [switch]$Refetch
)

$ErrorActionPreference = "Stop"
$sw = [Diagnostics.Stopwatch]::StartNew()

function Step($n, $msg) {
    Write-Host ""
    Write-Host "=== STEP $n : $msg ===" -ForegroundColor Cyan
}

function Need($path) {
    if (-not (Test-Path $path)) { throw "missing required input: $path" }
}

# ---- preflight -------------------------------------------------------------
Step 0 "preflight"
foreach ($f in @("model_scope_v1.txt", "decisions_final.csv",
                 "rules_v3_unfitted.txt", "build_network.py",
                 "make_network_v2.py", "make_rules.py",
                 "attractor_analysis.py", "run_screen.py")) { Need $f }

if (-not (Test-Path "run01\raw")) {
    Write-Warning "run01\raw not found -- OmniPath will be downloaded fresh."
    Write-Warning "The result may differ from the archived run, because"
    Write-Warning "OmniPath is updated continuously."
}

& $Python --version

# ---- 1. network from OmniPath ---------------------------------------------
Step 1 "build network from OmniPath  ->  final\"
$cache = if ($Refetch) { "final\raw" } else { "run01\raw" }
& $Python build_network.py `
    --anchors KRAS MDM2 `
    --sinks CASP3 CASP7 PARP1 CDKN1A `
    --nodes model_scope_v1.txt `
    --min-resources 2 --min-references 1 --min-curation-effort 5 `
    --terminal-sinks `
    --decisions decisions_final.csv `
    --cache $cache -d final
if ($LASTEXITCODE -ne 0) { throw "step 1 failed" }

# ---- 2. literature corrections --------------------------------------------
Step 2 "apply literature corrections  ->  network_v2_*.csv"
& $Python make_network_v2.py --in final\edges.csv --out .
if ($LASTEXITCODE -ne 0) { throw "step 2 failed" }

# ---- 3. validate the rules against the network ----------------------------
Step 3 "validate rules against network"
& $Python make_rules.py `
    --edges network_v2_edges.csv --nodes network_v2_nodes.csv `
    --check rules_v3_unfitted.txt `
    --fixed KRAS MDM2 MDM4 `
    --export-maboss model_v3.bnd --export-boolnet model_v3.bnet
if ($LASTEXITCODE -ne 0) {
    throw "step 3 failed -- rules disagree with the network; see output above"
}

# ---- 4. exact attractor analysis ------------------------------------------
Step 4 "exact attractor analysis  ->  attractor_report.md"
& $Python attractor_analysis.py `
    --rules rules_v3_unfitted.txt `
    --clamp control:KRAS=1,MDM2=1,MDM4=1 `
    --clamp MDM2_degraded:KRAS=1,MDM2=0,MDM4=1 `
    --clamp KRAS_degraded:KRAS=0,MDM2=1,MDM4=1 `
    --clamp both_degraded:KRAS=0,MDM2=0,MDM4=1 `
    --samples 120 --steps 6000 `
    --out attractor_report.md
if ($LASTEXITCODE -ne 0) { throw "step 4 failed" }

# ---- 5. rule-variant screen -----------------------------------------------
Step 5 "32-variant screen  ->  screen_report.md"
& $Python run_screen.py `
    --rules rules_v3_unfitted.txt `
    --out-json screen_results.json `
    --out-md screen_report.md
if ($LASTEXITCODE -ne 0) { throw "step 5 failed" }

# ---- checksums -------------------------------------------------------------
Step 6 "checksums of the report artifacts"
$artifacts = @("network_v2_edges.csv", "network_v2_nodes.csv",
               "network_v2_changelog.csv", "rules_v3_unfitted.txt",
               "attractor_report.md", "screen_report.md",
               "screen_results.json", "decisions_final.csv",
               "model_scope_v1.txt")
# NB: in Windows PowerShell 5.1 Get-FileHash does not accept bare strings
# from the pipeline (-Path binds by property name, not by value), so the
# files are passed explicitly.
$hashes = foreach ($f in $artifacts) {
    if (Test-Path $f) { Get-FileHash -Path $f -Algorithm SHA256 }
}
$hashes |
    Select-Object @{n = "file"; e = { Split-Path $_.Path -Leaf } }, Hash |
    Format-Table -AutoSize |
    Out-String |
    Tee-Object -FilePath checksums.txt |
    Write-Host

$sw.Stop()
Write-Host ""
Write-Host "PIPELINE COMPLETE in $([int]$sw.Elapsed.TotalSeconds)s" -ForegroundColor Green
Write-Host ""
Write-Host "Read in this order:"
Write-Host "  network_v2_changelog.csv   what changed from the OmniPath network and why"
Write-Host "  attractor_report.md        the four attractors, exact"
Write-Host "  screen_report.md           32-variant grid and predictions"
Write-Host ""
Write-Host "Note: line endings differ between Windows and Linux, so file hashes"
Write-Host "will not match the container's. Compare file CONTENT, not hashes."
