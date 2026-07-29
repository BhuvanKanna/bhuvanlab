#Requires -Version 5.1
<#
.SYNOPSIS
    Register (or remove) the auto-sync watcher as a Scheduled Task that starts at logon.

.DESCRIPTION
    Creates a task named "BhuvanlabAutoSync" that launches tools\sync-watch.ps1
    hidden, at every logon, and restarts it if it ever dies. The task runs in
    the current user's interactive session, so it inherits the GitHub
    credentials that `gh auth setup-git` configured -- no stored password needed.

.EXAMPLE
    .\tools\install-autosync.ps1
    Install and start it now.

.EXAMPLE
    .\tools\install-autosync.ps1 -Uninstall
    Stop and remove the scheduled task. The repo and scripts stay put.
#>
[CmdletBinding()]
param(
    # Remove the scheduled task instead of creating it.
    [switch]$Uninstall,

    # Quiet period handed to the watcher.
    [int]$QuietSeconds = 45,

    # Size guard handed to the watcher.
    [int]$ThresholdMB = 250
)

$ErrorActionPreference = 'Stop'

$TaskName  = 'BhuvanlabAutoSync'
$RepoRoot  = Split-Path -Parent $PSScriptRoot
$WatchPath = Join-Path $PSScriptRoot 'sync-watch.ps1'

# --- uninstall --------------------------------------------------------------
if ($Uninstall) {
    $existing = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    if ($null -eq $existing) {
        Write-Host "Task '$TaskName' is not registered; nothing to remove." -ForegroundColor Yellow
        exit 0
    }
    try { Stop-ScheduledTask -TaskName $TaskName -ErrorAction Stop } catch {}
    Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false
    Write-Host "Removed scheduled task '$TaskName'. Auto-sync is off." -ForegroundColor Green
    Write-Host "Manual syncing still works: .\tools\sync-now.ps1" -ForegroundColor Gray
    exit 0
}

# --- install ----------------------------------------------------------------
if (-not (Test-Path $WatchPath)) {
    Write-Host "Cannot find $WatchPath" -ForegroundColor Red
    exit 1
}

$argLine = '-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File "{0}" -QuietSeconds {1} -ThresholdMB {2}' -f `
           $WatchPath, $QuietSeconds, $ThresholdMB

$action = New-ScheduledTaskAction -Execute 'powershell.exe' `
                                  -Argument $argLine `
                                  -WorkingDirectory $RepoRoot

$trigger = New-ScheduledTaskTrigger -AtLogOn -User $env:USERNAME

# ExecutionTimeLimit 0 = never time out (it is a long-running watcher).
$settings = New-ScheduledTaskSettingsSet `
    -AllowStartIfOnBatteries `
    -DontStopIfGoingOnBatteries `
    -StartWhenAvailable `
    -RestartCount 3 `
    -RestartInterval (New-TimeSpan -Minutes 1) `
    -ExecutionTimeLimit (New-TimeSpan -Seconds 0)

$principal = New-ScheduledTaskPrincipal -UserId "$env:USERDOMAIN\$env:USERNAME" `
                                        -LogonType Interactive `
                                        -RunLevel Limited

$existing = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
if ($null -ne $existing) {
    Write-Host "Task '$TaskName' already exists -- replacing it." -ForegroundColor Yellow
    try { Stop-ScheduledTask -TaskName $TaskName -ErrorAction Stop } catch {}
    Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false
}

Register-ScheduledTask -TaskName $TaskName `
                       -Action $action `
                       -Trigger $trigger `
                       -Settings $settings `
                       -Principal $principal `
                       -Description 'Auto-commit and push bhuvanlab repo changes to GitHub.' | Out-Null

Start-ScheduledTask -TaskName $TaskName

Write-Host ""
Write-Host "Auto-sync installed and running." -ForegroundColor Green
Write-Host "  Task name     : $TaskName"
Write-Host "  Watching      : $RepoRoot"
Write-Host "  Quiet period  : ${QuietSeconds}s after the last write"
Write-Host "  Size guard    : ${ThresholdMB} MB (larger payloads commit locally and wait)"
Write-Host "  Log           : $(Join-Path $RepoRoot '.sync\sync.log')"
Write-Host ""
Write-Host "It restarts automatically at every logon." -ForegroundColor Gray
Write-Host "Turn it off with: .\tools\install-autosync.ps1 -Uninstall" -ForegroundColor Gray
