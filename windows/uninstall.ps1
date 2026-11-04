# Remove the LysoQuant server installed by install.ps1.
#
#   .\uninstall.ps1 [-Purge] [-AllUsers]     -Purge also deletes the process folder
#                                              -AllUsers removes the all-users installation
param([switch]$Purge, [switch]$AllUsers)
$ErrorActionPreference = "Stop"

if ($AllUsers) {
    $Task = "LysoQuant Server (all users)"
    $Prefix = Join-Path $env:ProgramFiles "LysoQuant"
    $Data = Join-Path $env:ProgramData "LysoQuant"
} else {
    $Task = "LysoQuant Server"
    $Prefix = Join-Path $env:LOCALAPPDATA "LysoQuant"
    $Data = $Prefix
}

# The process folder is the last argument of the server command
$Jobs = $null
$Cmd = Join-Path $Prefix "run-server.cmd"
if (Test-Path $Cmd) {
    $m = [regex]::Match((Get-Content $Cmd -Raw), '"([^"]*)" >> ')
    if ($m.Success) { $Jobs = $m.Groups[1].Value }
}

if (Get-ScheduledTask -TaskName $Task -ErrorAction SilentlyContinue) {
    Stop-ScheduledTask -TaskName $Task -ErrorAction SilentlyContinue
    Unregister-ScheduledTask -TaskName $Task -Confirm:$false
}
Get-Process -ErrorAction SilentlyContinue |
    Where-Object { $_.Path -and $_.Path.StartsWith($Prefix, [StringComparison]::OrdinalIgnoreCase) } |
    Stop-Process -Force
Start-Sleep -Seconds 2
Remove-Item $Prefix -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item (Join-Path $Data "LysoQuant.log") -Force -ErrorAction SilentlyContinue

if ($Purge -and $Jobs) { Remove-Item $Jobs -Recurse -Force -ErrorAction SilentlyContinue }
if ($AllUsers) { Remove-Item $Data -Force -ErrorAction SilentlyContinue }     # only if empty
Write-Host "LysoQuant server removed."
