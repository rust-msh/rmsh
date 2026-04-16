param(
    [switch]$SkipSmokeTest
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$targetDll = Join-Path $repoRoot "target\debug\_rmsh.dll"
$targetPyd = Join-Path $repoRoot "crates\py\python\rmsh\_rmsh.pyd"
$venvPython = Join-Path $repoRoot "crates\py\.venv\Scripts\python.exe"

Write-Host "[1/4] Building rmsh-py (debug) ..."
Push-Location $repoRoot
try {
    cargo build -p rmsh-py
} finally {
    Pop-Location
}

if (-not (Test-Path $targetDll)) {
    throw "Build finished but _rmsh.dll was not found at $targetDll"
}

Write-Host "[2/4] Refreshing Python extension ..."
Copy-Item -Force $targetDll $targetPyd

Write-Host "[3/4] Extension timestamp:"
Get-Item $targetPyd | Select-Object FullName, Length, LastWriteTime | Format-Table -AutoSize

if ($SkipSmokeTest) {
    Write-Host "[4/4] Smoke test skipped."
    exit 0
}

if (-not (Test-Path $venvPython)) {
    throw "Python venv not found at $venvPython"
}

Write-Host "[4/4] Running analytic primitive STEP smoke tests ..."
Push-Location $repoRoot
try {
    & $venvPython -c "import importlib, pathlib, rmsh; m=importlib.import_module('rmsh._rmsh'); print('using', m.__file__); rmsh.initialize(); rmsh.model.occ.addRectangle(0,0,0,2.0,1.0); rmsh.write('rect_runtime_check.step'); rmsh.clear(); rmsh.initialize(); rmsh.model.occ.addBox(0,0,0,1,2,3); rmsh.write('box_runtime_check.step'); rmsh.clear(); rmsh.initialize(); rmsh.model.occ.addSphere(0,0,0,1.0); rmsh.write('sphere_runtime_check.step'); rmsh.clear(); rmsh.initialize(); rmsh.model.occ.addCylinder(0,0,0,0,0,2,0.5); rmsh.write('cylinder_runtime_check.step'); rmsh.clear(); rmsh.initialize(); rmsh.model.occ.addCone(0,0,0,0,0,2,0.8,0.2); rmsh.write('cone_runtime_check.step'); rmsh.clear(); rmsh.initialize(); rmsh.model.occ.addTorus(0,0,0,0,0,1,1.0,0.2); rmsh.write('torus_runtime_check.step'); rmsh.finalize(); print(pathlib.Path('torus_runtime_check.step').resolve())"

    $files = @(
        "rect_runtime_check.step",
        "box_runtime_check.step",
        "sphere_runtime_check.step",
        "cylinder_runtime_check.step",
        "cone_runtime_check.step",
        "torus_runtime_check.step"
    )

    foreach ($f in $files) {
        if (-not (Test-Path $f)) {
            throw "Smoke test failed: $f was not generated"
        }

        $matches = Select-String -Path $f -Pattern "ADVANCED_FACE|PLANE|SPHERICAL_SURFACE|CYLINDRICAL_SURFACE|CONICAL_SURFACE|TOROIDAL_SURFACE|tri_face|tri_plane|TRIANGULATED|TESSELLATED" | ForEach-Object { $_.Line.Trim() }
        $joined = ($matches -join "`n")
        Write-Host "---- STEP markers: $f ----"
        Write-Host $joined

        if ($joined -match "tri_face|tri_plane|TRIANGULATED|TESSELLATED") {
            throw "Smoke test failed: fallback/tessellated markers still present in $f"
        }
        if ($joined -notmatch "ADVANCED_FACE") {
            throw "Smoke test failed: ADVANCED_FACE not found in $f"
        }
    }
} finally {
    Pop-Location
}

Write-Host "Done: rmsh Python extension refreshed and analytic primitive STEP smoke tests passed."
