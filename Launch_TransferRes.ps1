$ErrorActionPreference="Stop"
Set-Location $PSScriptRoot
uv run --with numpy --with pandas --with numba --with matplotlib python "$PSScriptRoot\TransferRes_GUI.py"
