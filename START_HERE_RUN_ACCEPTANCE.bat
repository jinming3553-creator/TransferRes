@echo off
cd /d "%~dp0"
echo TransferRes Windows Acceptance Hotfix v0.7.1
echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0Run_Acceptance_Test.ps1"
set "EXITCODE=%ERRORLEVEL%"
echo.
if "%EXITCODE%"=="0" (
  echo Acceptance finished: PASS
) else (
  echo Acceptance finished: FAIL
)
echo.
echo Please send ACCEPTANCE_REPORT.txt if you want me to review the result.
pause
exit /b %EXITCODE%
