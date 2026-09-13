$ErrorActionPreference="Stop"
Set-Location $PSScriptRoot

$version = "0.7.0-rc1"
$report = Join-Path $PSScriptRoot "RELEASE_BUILD_REPORT.txt"

"TransferRes Release Build" | Out-File $report -Encoding utf8
("Version: " + $version) | Out-File $report -Append -Encoding utf8
("Started: " + (Get-Date).ToString("s")) | Out-File $report -Append -Encoding utf8

Write-Host "Step 1/4: Running Windows acceptance..."
& powershell -ExecutionPolicy Bypass -File "$PSScriptRoot\Run_Acceptance_Test.ps1"
if ($LASTEXITCODE -ne 0) {
    "Acceptance: FAIL" | Out-File $report -Append -Encoding utf8
    throw "Acceptance failed. Release package will not be created."
}
"Acceptance: PASS" | Out-File $report -Append -Encoding utf8

$exe = Join-Path $PSScriptRoot "dist\TransferRes.exe"
if (-not (Test-Path $exe)) {
    throw "Expected executable missing: $exe"
}

Write-Host "Step 2/4: Preparing portable release folder..."
$portable = Join-Path $PSScriptRoot ("TransferRes_Portable_" + $version)
if (Test-Path $portable) { Remove-Item $portable -Recurse -Force }
New-Item -ItemType Directory -Force -Path $portable | Out-Null

Copy-Item $exe (Join-Path $portable "TransferRes.exe")
Copy-Item "$PSScriptRoot\QUICK_START.md" $portable
Copy-Item "$PSScriptRoot\MANUAL_GUI_ACCEPTANCE.md" $portable
Copy-Item "$PSScriptRoot\LICENSE" $portable
Copy-Item "$PSScriptRoot\CITATION.cff" $portable
Copy-Item "$PSScriptRoot\CHANGELOG.md" $portable
Copy-Item "$PSScriptRoot\REPRODUCIBILITY.md" $portable
Copy-Item "$PSScriptRoot\VERSION.txt" $portable
Copy-Item "$PSScriptRoot\ACCEPTANCE_REPORT.txt" $portable

Write-Host "Step 3/4: Generating checksum..."
$hash = (Get-FileHash (Join-Path $portable "TransferRes.exe") -Algorithm SHA256).Hash
$hash | Out-File (Join-Path $portable "TransferRes.exe.sha256.txt") -Encoding ascii
("SHA256: " + $hash) | Out-File $report -Append -Encoding utf8

Write-Host "Step 4/4: Creating portable ZIP..."
$zip = Join-Path $PSScriptRoot ("TransferRes_Portable_" + $version + ".zip")
if (Test-Path $zip) { Remove-Item $zip -Force }
Compress-Archive -Path (Join-Path $portable "*") -DestinationPath $zip -CompressionLevel Optimal

("Portable ZIP: " + $zip) | Out-File $report -Append -Encoding utf8
("Finished: " + (Get-Date).ToString("s")) | Out-File $report -Append -Encoding utf8
"OVERALL: PASS" | Out-File $report -Append -Encoding utf8

Write-Host ""
Write-Host "============================================================"
Write-Host "RELEASE BUILD PASS"
Write-Host "Portable ZIP:"
Write-Host $zip
Write-Host "============================================================"
