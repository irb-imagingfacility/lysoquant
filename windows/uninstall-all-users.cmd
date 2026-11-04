@echo off
rem Asks for administrator rights, then removes the all-users LysoQuant server
powershell -NoProfile -Command "Start-Process powershell -Verb RunAs -ArgumentList '-NoProfile -ExecutionPolicy Bypass -NoExit -File ""%~dp0uninstall.ps1"" -AllUsers %*'"
