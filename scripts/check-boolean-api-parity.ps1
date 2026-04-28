param(
    [switch]$CoreOnly,
    [switch]$StrictFragment,
    [switch]$SkipBuild
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$pythonExe = Join-Path $repoRoot "crates\py\.venv\Scripts\python.exe"
$compareScript = Join-Path $repoRoot "crates\py\examples\compare_boolean_api_contract_rmsh_gmsh.py"

if (-not (Test-Path $pythonExe)) {
    throw "Python venv not found at $pythonExe"
}

if (-not $SkipBuild) {
    Write-Host "[1/2] Refreshing rmsh Python extension ..."
    & powershell -ExecutionPolicy Bypass -File (Join-Path $repoRoot "scripts\refresh-rmsh-py.ps1") -SkipSmokeTest
}

Write-Host "[2/2] Running boolean API parity check ..."
$args = @($compareScript)
if ($CoreOnly) {
    $args += "--core-only"
}
if ($StrictFragment) {
    $args += "--strict-fragment"
}

& $pythonExe @args
