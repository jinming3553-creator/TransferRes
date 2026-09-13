@echo off
setlocal
cd /d "%~dp0"

where uv >nul 2>nul
if errorlevel 1 (
    echo.
    echo TransferRes needs "uv" for automatic launch.
    echo Install uv from https://docs.astral.sh/uv/
    echo.
    pause
    exit /b 1
)

echo Starting TransferRes...
uv run --with numpy --with pandas --with numba --with matplotlib python "%~dp0TransferRes_GUI.py"

if errorlevel 1 (
    echo.
    echo TransferRes exited with an error.
    pause
)
endlocal
