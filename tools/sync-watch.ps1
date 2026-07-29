#Requires -Version 5.1
<#
.SYNOPSIS
    Watch this repo and auto-commit + push whenever files change.

.DESCRIPTION
    Uses a FileSystemWatcher for instant detection, then waits for a quiet
    period (-QuietSeconds, default 45) before syncing. The quiet period matters:
    regenerating the fourparam tables writes for minutes, and committing
    mid-write would push a truncated CSV. Nothing is committed until the writes
    have actually stopped.

    Each sync delegates to sync-now.ps1, so the size guard and the data/ guard
    described there apply here too.

    Runs until stopped (Ctrl+C, or Stop-ScheduledTask). Logs to .sync\sync.log.

.EXAMPLE
    .\tools\sync-watch.ps1
    Watch with the default 45-second quiet period.

.EXAMPLE
    .\tools\sync-watch.ps1 -QuietSeconds 120
    Wait 2 minutes of silence before syncing -- better for long batch runs.
#>
[CmdletBinding()]
param(
    # Seconds of no filesystem activity required before a sync fires.
    [int]$QuietSeconds = 45,

    # Passed through to sync-now.ps1: hold pushes larger than this.
    [int]$ThresholdMB = 250,

    # Push large payloads and data/ changes without holding.
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot
$SyncNow  = Join-Path $PSScriptRoot 'sync-now.ps1'
$LogDir   = Join-Path $RepoRoot '.sync'
$LogFile  = Join-Path $LogDir 'sync.log'

if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir -Force | Out-Null }

function Write-Log {
    param([string]$Message, [string]$Level = 'INFO')
    $line = "{0}  [{1}]  {2}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $Level, $Message
    Add-Content -Path $LogFile -Value $line -Encoding utf8
    $color = 'Gray'
    if ($Level -eq 'OK')    { $color = 'Green' }
    if ($Level -eq 'WATCH') { $color = 'Cyan' }
    if ($Level -eq 'ERROR') { $color = 'Red' }
    Write-Host $line -ForegroundColor $color
}

# Directories whose contents must never trigger a sync. Matched per path
# SEGMENT, not as a substring: committing writes to .git and logging writes to
# .sync, and FileSystemWatcher reports those as directory-level events whose
# path is "...\.git" with no trailing separator. A substring test for "\.git\"
# misses that and the watcher then retriggers itself forever.
$IgnoredSegments = @('.git', '.sync', '__pycache__', '.venv', 'venv', 'env',
                     'fourparam-venv', '.vscode', '.idea')

$IgnoredExtensions = @('.pyc', '.pyo', '.tmp', '.swp', '.log')

function Test-Relevant {
    param([string]$FullPath)
    if ([string]::IsNullOrWhiteSpace($FullPath)) { return $false }

    $rel = $FullPath
    if ($FullPath.StartsWith($RepoRoot, [StringComparison]::OrdinalIgnoreCase)) {
        $rel = $FullPath.Substring($RepoRoot.Length)
    }
    $rel = $rel.Trim([char]'\', [char]'/')

    # The repo root directory itself -- fires whenever any subdirectory changes.
    if ($rel -eq '') { return $false }

    foreach ($segment in ($rel -split '[\\/]+')) {
        if ($IgnoredSegments -contains $segment) { return $false }
    }
    foreach ($ext in $IgnoredExtensions) {
        if ($rel.EndsWith($ext, [StringComparison]::OrdinalIgnoreCase)) { return $false }
    }
    if ($rel.EndsWith('~')) { return $false }

    return $true
}

# --- set up the watcher -----------------------------------------------------
$watcher = New-Object System.IO.FileSystemWatcher
$watcher.Path                  = $RepoRoot
$watcher.IncludeSubdirectories = $true
$watcher.NotifyFilter = [System.IO.NotifyFilters]::FileName `
                    -bor [System.IO.NotifyFilters]::DirectoryName `
                    -bor [System.IO.NotifyFilters]::LastWrite
$watcher.InternalBufferSize = 65536   # large regens generate a lot of events

$subs = @()
foreach ($eventName in @('Changed', 'Created', 'Deleted', 'Renamed')) {
    $subs += Register-ObjectEvent -InputObject $watcher -EventName $eventName `
                                  -SourceIdentifier ("RepoSync_" + $eventName)
}
$watcher.EnableRaisingEvents = $true

Write-Log "Watching $RepoRoot" 'WATCH'
Write-Log "Quiet period ${QuietSeconds}s | size guard ${ThresholdMB} MB | Ctrl+C to stop" 'WATCH'

# Reconcile once at startup. A FileSystemWatcher only reports events that happen
# while it is live, so anything edited before logon -- or while the watcher was
# stopped -- would otherwise sit uncommitted forever, invisible to the watcher.
function Invoke-Sync {
    param([string]$Reason)
    try {
        $syncArgs = @(
            '-NoProfile', '-ExecutionPolicy', 'Bypass', '-File', $SyncNow,
            '-ThresholdMB', $ThresholdMB, '-Quiet'
        )
        if ($Force) { $syncArgs += '-Force' }
        & powershell.exe @syncArgs
        # sync-now.ps1 logs its own outcome; 3 = intentionally held.
        if ($LASTEXITCODE -eq 3) {
            Write-Log "Sync held (see log above). Release: tools\sync-now.ps1 -Force" 'WATCH'
        }
    } catch {
        Write-Log ("Sync run failed (" + $Reason + "): " + $_.Exception.Message) 'ERROR'
    }
}

Write-Log "Startup reconciliation ..." 'WATCH'
Invoke-Sync 'startup'

$pending    = $false
$lastChange = Get-Date

try {
    while ($true) {

        # Drain every queued filesystem event.
        $sawRelevant = $false
        $queued = @(Get-Event -ErrorAction SilentlyContinue |
                    Where-Object { $_.SourceIdentifier -like 'RepoSync_*' })
        foreach ($evt in $queued) {
            $changedPath = $null
            if ($null -ne $evt.SourceEventArgs) { $changedPath = $evt.SourceEventArgs.FullPath }
            if (Test-Relevant $changedPath) { $sawRelevant = $true }
            Remove-Event -EventIdentifier $evt.EventIdentifier -ErrorAction SilentlyContinue
        }

        if ($sawRelevant) {
            if (-not $pending) { Write-Log "Change detected; waiting for writes to settle ..." }
            $pending    = $true
            $lastChange = Get-Date
        }

        # Fire once the writes have stopped for long enough.
        if ($pending -and ((Get-Date) - $lastChange).TotalSeconds -ge $QuietSeconds) {
            $pending = $false
            Invoke-Sync 'change'
        }

        Start-Sleep -Seconds 3
    }
}
finally {
    $watcher.EnableRaisingEvents = $false
    foreach ($sub in $subs) {
        Unregister-Event -SubscriptionId $sub.Id -ErrorAction SilentlyContinue
    }
    $watcher.Dispose()
    Write-Log "Watcher stopped." 'WATCH'
}
