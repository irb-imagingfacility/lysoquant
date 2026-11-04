@echo off
rem Asks for administrator rights, then installs the LysoQuant server for all users of this PC
powershell -NoProfile -Command "Start-Process powershell -Verb RunAs -ArgumentList '-NoProfile -ExecutionPolicy Bypass -NoExit -File ""%~dp0install.ps1"" -AllUsers %*'"
