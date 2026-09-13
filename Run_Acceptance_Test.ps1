$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$report = Join-Path $PSScriptRoot "ACCEPTANCE_REPORT.txt"
$log = Join-Path $PSScriptRoot "ACCEPTANCE_BUILD_LOG.txt"
$selfLog = Join-Path $PSScriptRoot "ACCEPTANCE_SELF_TEST_LOG.txt"

"TransferRes Windows Acceptance Test" | Out-File $report -Encoding utf8
("Started: " + (Get-Date).ToString("s")) | Out-File $report -Append -Encoding utf8
"" | Out-File $report -Append -Encoding utf8

function Add-Result($name, $status, $detail="") {
    $line = "[$status] $name"
    if ($detail -ne "") { $line += " - $detail" }
    $line | Out-File $report -Append -Encoding utf8
    Write-Host $line
}

function Invoke-ProcessCaptured {
    param(
        [Parameter(Mandatory=$true)][string]$FilePath,
        [Parameter(Mandatory=$false)][string[]]$ArgumentList = @(),
        [Parameter(Mandatory=$true)][string]$StdoutPath,
        [Parameter(Mandatory=$true)][string]$StderrPath
    )

    if (Test-Path $StdoutPath) { Remove-Item $StdoutPath -Force }
    if (Test-Path $StderrPath) { Remove-Item $StderrPath -Force }

    $p = Start-Process `
        -FilePath $FilePath `
        -ArgumentList $ArgumentList `
        -WorkingDirectory $PSScriptRoot `
        -RedirectStandardOutput $StdoutPath `
        -RedirectStandardError $StderrPath `
        -NoNewWindow `
        -Wait `
        -PassThru

    return $p.ExitCode
}

try {
    # 1. Environment
    if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
        Add-Result "uv available" "FAIL" "uv not found on PATH"
        throw "uv is required."
    }

    $uvVersion = (& uv --version 2>$null | Out-String).Trim()
    if ($LASTEXITCODE -ne 0) {
        Add-Result "uv available" "FAIL" "uv --version failed"
        throw "uv is not callable."
    }
    Add-Result "uv available" "PASS" $uvVersion

    # 2. Build
    Write-Host ""
    Write-Host "Building standalone EXE..."
    Write-Host "(Download/progress text on stderr is normal and is NOT treated as failure.)"

    $buildOut = Join-Path $PSScriptRoot "_accept_build_stdout.tmp"
    $buildErr = Join-Path $PSScriptRoot "_accept_build_stderr.tmp"

    $buildCode = Invoke-ProcessCaptured `
        -FilePath "powershell.exe" `
        -ArgumentList @(
            "-NoProfile",
            "-ExecutionPolicy", "Bypass",
            "-File", (Join-Path $PSScriptRoot "Build_TransferRes_EXE.ps1")
        ) `
        -StdoutPath $buildOut `
        -StderrPath $buildErr

    # Combine logs after process completion. Do not pipe native stderr into a
    # Stop-on-error PowerShell pipeline.
    "=== STDOUT ===" | Out-File $log -Encoding utf8
    if (Test-Path $buildOut) {
        Get-Content $buildOut | Out-File $log -Append -Encoding utf8
        Get-Content $buildOut | ForEach-Object { Write-Host $_ }
    }
    "=== STDERR (may contain normal uv progress) ===" | Out-File $log -Append -Encoding utf8
    if (Test-Path $buildErr) {
        Get-Content $buildErr | Out-File $log -Append -Encoding utf8
        Get-Content $buildErr | ForEach-Object { Write-Host $_ }
    }

    Remove-Item $buildOut,$buildErr -Force -ErrorAction SilentlyContinue

    if ($buildCode -ne 0) {
        Add-Result "PyInstaller build" "FAIL" "Build returned exit code $buildCode"
        throw "Build failed. See ACCEPTANCE_BUILD_LOG.txt."
    }
    Add-Result "PyInstaller build" "PASS"

    $exe = Join-Path $PSScriptRoot "dist\TransferRes.exe"
    if (-not (Test-Path $exe)) {
        Add-Result "EXE created" "FAIL" "dist\TransferRes.exe missing"
        throw "Executable missing."
    }
    Add-Result "EXE created" "PASS" $exe

    $sizeMB = [math]::Round((Get-Item $exe).Length / 1MB, 2)
    Add-Result "EXE size readable" "PASS" "$sizeMB MB"

    # 3. Packaged EXE self-test
    Write-Host ""
    Write-Host "Running packaged EXE self-test..."

    $selfOutFile = Join-Path $PSScriptRoot "_accept_self_stdout.tmp"
    $selfErrFile = Join-Path $PSScriptRoot "_accept_self_stderr.tmp"

    $selfCode = Invoke-ProcessCaptured `
        -FilePath $exe `
        -ArgumentList @("--self-test") `
        -StdoutPath $selfOutFile `
        -StderrPath $selfErrFile

    "=== STDOUT ===" | Out-File $selfLog -Encoding utf8
    if (Test-Path $selfOutFile) {
        Get-Content $selfOutFile | Out-File $selfLog -Append -Encoding utf8
    }
    "=== STDERR ===" | Out-File $selfLog -Append -Encoding utf8
    if (Test-Path $selfErrFile) {
        Get-Content $selfErrFile | Out-File $selfLog -Append -Encoding utf8
    }

    $selfOut = ""
    if (Test-Path $selfOutFile) {
        $selfOut = Get-Content $selfOutFile -Raw
    }
    if (Test-Path $selfErrFile) {
        $selfOut += "`n" + (Get-Content $selfErrFile -Raw)
    }

    Remove-Item $selfOutFile,$selfErrFile -Force -ErrorAction SilentlyContinue

    if ($selfCode -ne 0) {
        Add-Result "Packaged EXE self-test" "FAIL" "Exit code $selfCode"
        throw "Packaged executable self-test failed."
    }

    if ($selfOut -notmatch "TRANSFERRES_SELF_TEST=PASS") {
        Add-Result "Packaged EXE self-test marker" "FAIL" "PASS marker not found"
        throw "PASS marker missing."
    }
    Add-Result "Packaged EXE self-test" "PASS"

    # 4. Frozen regressions
    if ($selfOut -match "SELECTED_LEVEL=class_name") {
        Add-Result "Frozen sample selected level" "PASS" "class_name"
    } else {
        Add-Result "Frozen sample selected level" "FAIL"
        throw "Selected-level regression failed."
    }

    if ($selfOut -match "ADJUSTED_FEATURE_COUNTS=\[20\]") {
        Add-Result "Feature-count overflow handling" "PASS" "500/1000/2000 -> 20"
    } else {
        Add-Result "Feature-count overflow handling" "FAIL"
        throw "Feature-count handling regression failed."
    }

    # 5. SHA256
    $hash = (Get-FileHash $exe -Algorithm SHA256).Hash
    Add-Result "SHA256 generated" "PASS" $hash
    $hash | Out-File (Join-Path $PSScriptRoot "TransferRes.exe.sha256.txt") -Encoding ascii

    # 6. Path with spaces
    $spaceDir = Join-Path $env:TEMP "TransferRes Acceptance Test"
    New-Item -ItemType Directory -Force -Path $spaceDir | Out-Null
    $spaceExe = Join-Path $spaceDir "TransferRes.exe"
    Copy-Item $exe $spaceExe -Force

    $spOut = Join-Path $spaceDir "stdout.txt"
    $spErr = Join-Path $spaceDir "stderr.txt"
    $spaceCode = Invoke-ProcessCaptured `
        -FilePath $spaceExe `
        -ArgumentList @("--self-test") `
        -StdoutPath $spOut `
        -StderrPath $spErr

    $spaceText = ""
    if (Test-Path $spOut) { $spaceText += (Get-Content $spOut -Raw) }
    if (Test-Path $spErr) { $spaceText += "`n" + (Get-Content $spErr -Raw) }

    if ($spaceCode -eq 0 -and $spaceText -match "TRANSFERRES_SELF_TEST=PASS") {
        Add-Result "Path containing spaces" "PASS"
    } else {
        Add-Result "Path containing spaces" "FAIL" "Exit code $spaceCode"
        throw "Space-path self-test failed."
    }

    # 7. Non-ASCII path
    $unicodeDir = Join-Path $env:TEMP "TransferRes_验收"
    New-Item -ItemType Directory -Force -Path $unicodeDir | Out-Null
    $unicodeExe = Join-Path $unicodeDir "TransferRes.exe"
    Copy-Item $exe $unicodeExe -Force

    $uniOut = Join-Path $unicodeDir "stdout.txt"
    $uniErr = Join-Path $unicodeDir "stderr.txt"
    $unicodeCode = Invoke-ProcessCaptured `
        -FilePath $unicodeExe `
        -ArgumentList @("--self-test") `
        -StdoutPath $uniOut `
        -StderrPath $uniErr

    $unicodeText = ""
    if (Test-Path $uniOut) { $unicodeText += (Get-Content $uniOut -Raw) }
    if (Test-Path $uniErr) { $unicodeText += "`n" + (Get-Content $uniErr -Raw) }

    if ($unicodeCode -eq 0 -and $unicodeText -match "TRANSFERRES_SELF_TEST=PASS") {
        Add-Result "Non-ASCII path" "PASS"
    } else {
        Add-Result "Non-ASCII path" "FAIL" "Exit code $unicodeCode"
        throw "Unicode-path self-test failed."
    }

    "" | Out-File $report -Append -Encoding utf8
    "OVERALL: PASS" | Out-File $report -Append -Encoding utf8
    ("Finished: " + (Get-Date).ToString("s")) | Out-File $report -Append -Encoding utf8

    Write-Host ""
    Write-Host "============================================================"
    Write-Host "OVERALL: PASS"
    Write-Host "Acceptance report:"
    Write-Host $report
    Write-Host "============================================================"
    exit 0
}
catch {
    "" | Out-File $report -Append -Encoding utf8
    "OVERALL: FAIL" | Out-File $report -Append -Encoding utf8
    ("Error: " + $_.Exception.Message) | Out-File $report -Append -Encoding utf8
    ("Finished: " + (Get-Date).ToString("s")) | Out-File $report -Append -Encoding utf8

    Write-Host ""
    Write-Host "============================================================"
    Write-Host "OVERALL: FAIL"
    Write-Host $_.Exception.Message
    Write-Host "Send me:"
    Write-Host $report
    Write-Host $log
    if (Test-Path $selfLog) { Write-Host $selfLog }
    Write-Host "============================================================"
    exit 1
}
