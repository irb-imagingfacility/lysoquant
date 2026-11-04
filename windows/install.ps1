# Install the LysoQuant server for the current user on Windows. It runs as a
# background task that starts when you log in; no administrator rights and no
# Python needed (a standalone Python is downloaded with uv).
#
#   .\install.ps1 [-Jobs DIR] [-Gpu "GPU 0"|none] [-AllUsers]
#
# or double-click install.cmd. -AllUsers (install-all-users.cmd) installs one
# server for every user of the PC instead: it needs administrator rights, runs
# as SYSTEM from startup (no login needed) and shares one process folder. Works from the repository (windows\install.ps1,
# uses ..\server and ..\model) and from the package made by build-package.sh
# (install.ps1 next to server\ and model\).
param(
    [string]$Jobs = "",          # default %USERPROFILE%\LysoQuant\jobs, or %ProgramData%\LysoQuant\jobs with -AllUsers
    [string]$Gpu = "GPU 0",      # "GPU 0" selects the NVIDIA GPU, "none" the CPU
    [switch]$AllUsers
)
$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"      # much faster downloads in Windows PowerShell
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

if ($AllUsers) {
    $IsAdmin = ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole(
        [Security.Principal.WindowsBuiltInRole]::Administrator)
    if (-not $IsAdmin) { throw "-AllUsers needs administrator rights: use install-all-users.cmd" }
    $Task = "LysoQuant Server (all users)"
    $Prefix = Join-Path $env:ProgramFiles "LysoQuant"    # only administrators can write here: the server runs as SYSTEM
    $Data = Join-Path $env:ProgramData "LysoQuant"
    if (-not $Jobs) { $Jobs = Join-Path $Data "jobs" }
} else {
    $Task = "LysoQuant Server"
    $Prefix = Join-Path $env:LOCALAPPDATA "LysoQuant"
    $Data = $Prefix
    if (-not $Jobs) { $Jobs = Join-Path $env:USERPROFILE "LysoQuant\jobs" }
}
$Log = Join-Path $Data "LysoQuant.log"
$PythonVersion = "3.12"
$CudaIndex = "https://download.pytorch.org/whl/cu124"     # PyPI only has CPU PyTorch for Windows

$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
if (Test-Path (Join-Path $Here "server")) { $Root = $Here } else { $Root = Split-Path -Parent $Here }
$ServerDir = Join-Path $Root "server"
$Model = Join-Path $Root "model\lysoquant-7-16.pt"
if (-not (Test-Path $Model)) { throw "Model $Model not found" }
if (-not (Test-Path $ServerDir)) { throw "Server sources $ServerDir not found" }

$arch = if ($env:PROCESSOR_ARCHITEW6432) { $env:PROCESSOR_ARCHITEW6432 } else { $env:PROCESSOR_ARCHITECTURE }
switch ($arch) {
    "AMD64" { $UvArch = "x86_64" }
    "ARM64" { $UvArch = "aarch64" }
    default { throw "Unsupported architecture $arch" }
}
$UseCuda = ($arch -eq "AMD64") -and ($Gpu -ne "none") -and [bool](Get-Command nvidia-smi -ErrorAction SilentlyContinue)

Write-Host "Installing into $Prefix"
New-Item -ItemType Directory -Force (Join-Path $Prefix "bin"), (Join-Path $Prefix "models"), $Data | Out-Null

# Stop a previous installation: Windows cannot replace files that are in use
if (Get-ScheduledTask -TaskName $Task -ErrorAction SilentlyContinue) {
    Stop-ScheduledTask -TaskName $Task -ErrorAction SilentlyContinue
}
Get-Process -ErrorAction SilentlyContinue |
    Where-Object { $_.Path -and $_.Path.StartsWith($Prefix, [StringComparison]::OrdinalIgnoreCase) } |
    Stop-Process -Force
Start-Sleep -Seconds 2

$Uv = Join-Path $Prefix "bin\uv.exe"
if (-not (Test-Path $Uv)) {
    Write-Host "Downloading uv"
    $tmp = Join-Path ([IO.Path]::GetTempPath()) ([IO.Path]::GetRandomFileName())
    New-Item -ItemType Directory $tmp | Out-Null
    try {
        $zip = Join-Path $tmp "uv.zip"
        Invoke-WebRequest "https://github.com/astral-sh/uv/releases/latest/download/uv-$UvArch-pc-windows-msvc.zip" -OutFile $zip
        Expand-Archive $zip -DestinationPath $tmp
        Copy-Item (Get-ChildItem $tmp -Recurse -Filter uv.exe | Select-Object -First 1).FullName $Uv
    } finally {
        Remove-Item $tmp -Recurse -Force -ErrorAction SilentlyContinue
    }
}

function Run($exe, $arguments) {
    & $exe @arguments
    if ($LASTEXITCODE -ne 0) { throw "$exe failed (exit code $LASTEXITCODE)" }
}

$env:UV_PYTHON_INSTALL_DIR = Join-Path $Prefix "python"
$env:UV_CACHE_DIR = Join-Path $Prefix "cache"
$Venv = Join-Path $Prefix "venv"
$Python = Join-Path $Venv "Scripts\python.exe"
Write-Host "Creating the Python environment (downloads PyTorch, a few hundred MB to 2.5 GB with CUDA)"
Run $Uv @("venv", "--python", $PythonVersion, "--clear", $Venv)
if ($UseCuda) {
    Write-Host "NVIDIA GPU found: installing PyTorch with CUDA"
    Run $Uv @("pip", "install", "--python", $Python, "torch", "--index-url", $CudaIndex, "--extra-index-url", "https://pypi.org/simple")
} else {
    Write-Host "No NVIDIA GPU selected: installing the CPU version of PyTorch"
}
Run $Uv @("pip", "install", "--python", $Python, $ServerDir)

Copy-Item $Model (Join-Path $Prefix "models\lysoquant-7-16.pt") -Force
New-Item -ItemType Directory -Force $Jobs | Out-Null
if ($AllUsers) {
    # Every local user can create, read and delete job folders (the server, as SYSTEM, can always read them)
    Run icacls @($Jobs, "/grant", "*S-1-5-32-545:(OI)(CI)M")
}

# Launcher: the task runs a hidden wscript, which runs the server with its output in the log
$Cmd = Join-Path $Prefix "run-server.cmd"
$Vbs = Join-Path $Prefix "run-server.vbs"
@"
@echo off
set PYTHONUNBUFFERED=1
"$Venv\Scripts\lysoquant-serve.exe" --weights "$Prefix\models\lysoquant-7-16.pt" --gpu "$Gpu" "$Jobs" >> "$Log" 2>&1
"@ | Set-Content $Cmd -Encoding ASCII
@"
Set sh = CreateObject("WScript.Shell")
WScript.Quit sh.Run("cmd /c ""$Cmd""", 0, True)
"@ | Set-Content $Vbs -Encoding ASCII

if ($AllUsers) {
    $action = New-ScheduledTaskAction -Execute $Cmd      # SYSTEM has no desktop: no window to hide
    $trigger = New-ScheduledTaskTrigger -AtStartup
    $principal = New-ScheduledTaskPrincipal -UserId "SYSTEM" -LogonType ServiceAccount -RunLevel Highest
} else {
    $action = New-ScheduledTaskAction -Execute "wscript.exe" -Argument "`"$Vbs`""
    $trigger = New-ScheduledTaskTrigger -AtLogOn -User "$env:USERDOMAIN\$env:USERNAME"
    $principal = New-ScheduledTaskPrincipal -UserId "$env:USERDOMAIN\$env:USERNAME" -LogonType Interactive
}
$settings = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
    -StartWhenAvailable -MultipleInstances IgnoreNew -ExecutionTimeLimit ([TimeSpan]::Zero) `
    -RestartCount 999 -RestartInterval (New-TimeSpan -Minutes 1)
Register-ScheduledTask -TaskName $Task -Action $action -Trigger $trigger -Settings $settings -Principal $principal `
    -Description "LysoQuant segmentation server" -Force | Out-Null
Start-ScheduledTask -TaskName $Task

$device = & $Python -c "import torch; print('CUDA: ' + torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU only (no CUDA)')"
Write-Host @"

LysoQuant server installed and running (starts at every $(if ($AllUsers) { "boot, for all users" } else { "login" })). PyTorch: $device
  Process folder  $Jobs   <- set it as "Process folder on server" in Fiji,
                             with "Use remote host" unticked
  Log             Get-Content "$Log" -Wait   ("Watching ... for jobs" once the model is loaded)
  Uninstall       $(Join-Path $Here $(if ($AllUsers) { "uninstall-all-users.cmd" } else { "uninstall.cmd" }))
"@
