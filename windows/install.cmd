@echo off
rem Double-click to install the LysoQuant server (no administrator rights needed)
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" %*
echo.
pause
