$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host "Building standalone Windows TransferRes.exe..."
Write-Host "This must be run on Windows because PyInstaller is not a cross-compiler."

# Important:
# uv/PyInstaller may write ordinary download/progress information to STDERR.
# PowerShell 5.1 can surface native STDERR as ErrorRecord objects when the
# caller redirects/merges streams. We therefore judge native success ONLY
# by the process exit code, not by whether STDERR contains text.

$uvArgs = @(
    "run",
    "--with", "pyinstaller",
    "--with", "numpy",
    "--with", "pandas",
    "--with", "numba",
    "--with", "matplotlib",
    "pyinstaller",
    "--noconfirm",
    "--clean",
    "--windowed",
    "--onefile",
    "--name", "TransferRes",
    "--hidden-import", "numba",
    "--hidden-import", "matplotlib.backends.backend_agg",
    "--add-data", "transferres;transferres",
    "--add-data", "sample_data;sample_data",
    "--add-data", "interpretation.py;.",
    "TransferRes_GUI.py"
)

& uv @uvArgs
$code = $LASTEXITCODE

if ($code -ne 0) {
    throw "uv/PyInstaller exited with code $code."
}

$exe = Join-Path $PSScriptRoot "dist\TransferRes.exe"
if (-not (Test-Path $exe)) {
    throw "PyInstaller returned success but dist\TransferRes.exe was not created."
}

Write-Host ""
Write-Host "DONE."
Write-Host "Standalone executable:"
Write-Host $exe
exit 0
