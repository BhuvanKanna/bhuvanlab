#Requires -Version 5.1
<#
.SYNOPSIS
    Commit and push any local changes in this repo to GitHub.

.DESCRIPTION
    Stages every change (git add -A stages ONLY altered and new files -- unchanged
    files are never re-uploaded), commits with a generated message, and pushes.

    Two situations cause a push to be HELD rather than sent automatically:

      1. The pending upload exceeds -ThresholdMB (default 250 MB raw). A full
         100-table regeneration is ~1.9 GB raw and would tie up the connection
         for a long time, so it waits for a deliberate go-ahead.
      2. Files under data/ changed. Those are Git LFS pointers and the repo is
         configured --skip-smudge, so unexpected changes there deserve a human
         look before they reach the remote.

    In both cases the commit is still made locally -- nothing is lost. Release it
    with:  .\tools\sync-now.ps1 -Force

.EXAMPLE
    .\tools\sync-now.ps1
    Normal sync. Pushes if under the size guard.

.EXAMPLE
    .\tools\sync-now.ps1 -Force
    Push regardless of size, and release any previously held commits.
#>
[CmdletBinding()]
param(
    # Push even when the size guard or the data/ guard would otherwise hold it.
    [switch]$Force,

    # Hold pushes whose raw payload exceeds this many megabytes.
    [int]$ThresholdMB = 250,

    # Suppress console output (used by the watcher; the log file is still written).
    [switch]$Quiet
)

$ErrorActionPreference = 'Stop'

$RepoRoot = Split-Path -Parent $PSScriptRoot
$LogDir   = Join-Path $RepoRoot '.sync'
$LogFile  = Join-Path $LogDir 'sync.log'

if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir -Force | Out-Null }

function Write-Log {
    param([string]$Message, [string]$Level = 'INFO')
    $line = "{0}  [{1}]  {2}" -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'), $Level, $Message
    Add-Content -Path $LogFile -Value $line -Encoding utf8
    if (-not $Quiet) {
        $color = 'Gray'
        if ($Level -eq 'OK')    { $color = 'Green' }
        if ($Level -eq 'HOLD')  { $color = 'Yellow' }
        if ($Level -eq 'ERROR') { $color = 'Red' }
        Write-Host $line -ForegroundColor $color
    }
}

# --- locate git (Task Scheduler may not have it on PATH) --------------------
$GitCmd = $null
$gitLookup = Get-Command git -ErrorAction SilentlyContinue
if ($null -ne $gitLookup) {
    $GitCmd = $gitLookup.Source
} else {
    foreach ($candidate in @(
        "$env:ProgramFiles\Git\cmd\git.exe",
        "${env:ProgramFiles(x86)}\Git\cmd\git.exe",
        "$env:LOCALAPPDATA\Programs\Git\cmd\git.exe"
    )) {
        if (Test-Path $candidate) { $GitCmd = $candidate; break }
    }
}
if ($null -eq $GitCmd) {
    Write-Log "git executable not found. Install Git or add it to PATH." 'ERROR'
    exit 1
}

# Invoke git in the repo, returning stdout lines. Never throws on non-zero.
function Invoke-Git {
    param([Parameter(ValueFromRemainingArguments = $true)][string[]]$GitArgs)
    $out = & $GitCmd -C $RepoRoot @GitArgs 2>&1
    $script:GitExit = $LASTEXITCODE
    return $out
}

# ---------------------------------------------------------------------------
# 1. Is there anything to do?
# ---------------------------------------------------------------------------
$porcelain = @(Invoke-Git status --porcelain)
if ($GitExit -ne 0) {
    Write-Log ("git status failed: " + ($porcelain -join ' ')) 'ERROR'
    exit 1
}

# ---------------------------------------------------------------------------
# 2. Stage and commit whatever changed
# ---------------------------------------------------------------------------
if ($porcelain.Count -gt 0) {

    Invoke-Git add -A | Out-Null
    if ($GitExit -ne 0) { Write-Log "git add failed" 'ERROR'; exit 1 }

    # Recount from the index -- .gitignore may have excluded everything we saw.
    $staged = @(Invoke-Git diff --cached --name-status)

    if ($staged.Count -eq 0) {
        Write-Log "Changes seen but all ignored by .gitignore; nothing to commit."
    } else {
        $added    = @($staged | Where-Object { $_ -match '^A' }).Count
        $modified = @($staged | Where-Object { $_ -match '^M' }).Count
        $deleted  = @($staged | Where-Object { $_ -match '^D' }).Count
        $renamed  = @($staged | Where-Object { $_ -match '^R' }).Count

        $parts = @()
        if ($added -gt 0)    { $parts += "$added added" }
        if ($modified -gt 0) { $parts += "$modified modified" }
        if ($deleted -gt 0)  { $parts += "$deleted deleted" }
        if ($renamed -gt 0)  { $parts += "$renamed renamed" }
        $summary = $parts -join ', '

        # Body: the changed paths (capped so a full regen doesn't make a 100-line message)
        $paths = $staged | ForEach-Object { ($_ -split "`t")[-1] }
        $body  = ($paths | Select-Object -First 20) -join "`n"
        if ($paths.Count -gt 20) {
            $body = $body + "`n... and " + ($paths.Count - 20) + " more"
        }

        $msg = "auto-sync: $summary`n`n$body"
        $commitOut = Invoke-Git commit -m $msg
        if ($GitExit -ne 0) {
            Write-Log ("git commit failed: " + ($commitOut -join ' ')) 'ERROR'
            exit 1
        }
        Write-Log "Committed: $summary"
    }
} else {
    Write-Log "No file changes."
}

# ---------------------------------------------------------------------------
# 3. Anything unpushed? (covers both a fresh commit and a previously held one)
# ---------------------------------------------------------------------------
Invoke-Git fetch --filter=blob:none origin main --quiet | Out-Null

$ahead = @(Invoke-Git rev-list --count 'origin/main..main')
$behind = @(Invoke-Git rev-list --count 'main..origin/main')
$aheadN  = 0
$behindN = 0
if ($ahead.Count -gt 0)  { [int]::TryParse($ahead[0].Trim(),  [ref]$aheadN)  | Out-Null }
if ($behind.Count -gt 0) { [int]::TryParse($behind[0].Trim(), [ref]$behindN) | Out-Null }

if ($aheadN -eq 0) {
    if ($behindN -gt 0) {
        Write-Log "Nothing to push. Remote is $behindN commit(s) ahead -- run: git pull --rebase" 'HOLD'
    } else {
        Write-Log "Already in sync with GitHub." 'OK'
    }
    exit 0
}

if ($behindN -gt 0) {
    Write-Log "Diverged: $aheadN local / $behindN remote commit(s). Not pushing automatically." 'HOLD'
    Write-Log "Resolve manually, e.g.:  git -C `"$RepoRoot`" pull --rebase" 'HOLD'
    exit 2
}

# ---------------------------------------------------------------------------
# 4. Size + data/ guards
# ---------------------------------------------------------------------------
$changedPaths = @(Invoke-Git diff --name-only 'origin/main..main')

$rawBytes = 0
foreach ($rel in $changedPaths) {
    if ([string]::IsNullOrWhiteSpace($rel)) { continue }
    $full = Join-Path $RepoRoot $rel
    if (Test-Path $full -PathType Leaf) {
        $rawBytes += (Get-Item $full).Length
    }
}
$rawMB = [math]::Round($rawBytes / 1MB, 1)

$dataTouched = @($changedPaths | Where-Object { $_ -like 'fourparamsacrosstissues/data/*' })

$holdReasons = @()
if ($rawMB -gt $ThresholdMB) {
    $holdReasons += "payload ${rawMB} MB exceeds ${ThresholdMB} MB guard"
}
if ($dataTouched.Count -gt 0) {
    $holdReasons += "$($dataTouched.Count) LFS pointer file(s) under data/ changed"
}

if ($holdReasons.Count -gt 0 -and -not $Force) {
    Write-Log ("Push HELD -- " + ($holdReasons -join '; ')) 'HOLD'
    Write-Log "$aheadN commit(s) are committed locally and safe. Release with:" 'HOLD'
    Write-Log "  powershell -File `"$PSCommandPath`" -Force" 'HOLD'
    exit 3
}

# ---------------------------------------------------------------------------
# 5. Push
# ---------------------------------------------------------------------------
Write-Log "Pushing $aheadN commit(s), ~${rawMB} MB raw ..."
$pushOut = Invoke-Git push origin main
if ($GitExit -ne 0) {
    Write-Log ("git push failed: " + ($pushOut -join ' ')) 'ERROR'
    exit 1
}
Write-Log "Pushed to github.com/BhuvanKanna/bhuvanlab (main). ~${rawMB} MB raw." 'OK'
exit 0
