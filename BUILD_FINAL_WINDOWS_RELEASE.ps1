$ErrorActionPreference="Stop"
Set-Location $PSScriptRoot
$version="1.0.0"

Write-Host "TransferRes $version final Windows release build"
Write-Host ""

& powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$PSScriptRoot\Run_Acceptance_Test.ps1"
if ($LASTEXITCODE -ne 0) {
    throw "Final Windows acceptance failed. Release ZIP was not created."
}

$exe=Join-Path $PSScriptRoot "dist\TransferRes.exe"
if (-not (Test-Path $exe)) { throw "dist\TransferRes.exe missing." }

$portable=Join-Path $PSScriptRoot ("TransferRes_" + $version + "_Windows_Portable")
if (Test-Path $portable) { Remove-Item $portable -Recurse -Force }
New-Item -ItemType Directory -Path $portable | Out-Null

$files=@(
  "README.md",
  "VERSION.txt",
  "CHANGELOG.md",
  "RELEASE_NOTES_v1.0.0.md",
  "REPRODUCIBILITY.md",
  "CITATION.cff",
  "LICENSE",
  "ACCEPTANCE_REPORT.txt"
)

Copy-Item $exe (Join-Path $portable "TransferRes.exe")
foreach($f in $files) {
    $p=Join-Path $PSScriptRoot $f
    if (Test-Path $p) { Copy-Item $p $portable }
}

$hash=(Get-FileHash (Join-Path $portable "TransferRes.exe") -Algorithm SHA256).Hash
$hash | Out-File (Join-Path $portable "TransferRes.exe.sha256.txt") -Encoding ascii

$zip=Join-Path $PSScriptRoot ("TransferRes_" + $version + "_Windows_Portable.zip")
if (Test-Path $zip) { Remove-Item $zip -Force }
Compress-Archive -Path (Join-Path $portable "*") -DestinationPath $zip -CompressionLevel Optimal

Write-Host ""
Write-Host "FINAL RELEASE BUILD: PASS"
Write-Host $zip
Write-Host ("EXE SHA256: " + $hash)
