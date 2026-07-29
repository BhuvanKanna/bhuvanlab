# Auto-sync to GitHub

Keeps this working copy continuously synced to
[`BhuvanKanna/bhuvanlab`](https://github.com/BhuvanKanna/bhuvanlab) on `main`.

Edit a file, wait ~45 seconds, and it is committed and pushed. Only altered and
new files are sent — unchanged files are never re-uploaded (that is native `git`
behaviour, not something these scripts have to arrange).

## Setup (already done on this machine)

```powershell
gh auth setup-git                      # credential helper, so pushes never prompt
.\tools\install-autosync.ps1           # register + start the watcher
```

The watcher is a Scheduled Task named `BhuvanlabAutoSync` that restarts at every
logon.

## The three scripts

| Script | What it does |
|---|---|
| `sync-now.ps1` | One-shot: stage → commit → push. Use it to sync immediately or to release a held push. |
| `sync-watch.ps1` | Long-running watcher. Detects changes, waits for quiet, calls `sync-now.ps1`. |
| `install-autosync.ps1` | Registers/removes the watcher as a Scheduled Task. |

## Everyday commands

```powershell
.\tools\sync-now.ps1                      # sync right now
.\tools\sync-now.ps1 -Force               # sync and release a held push
.\tools\install-autosync.ps1 -Uninstall   # turn auto-sync off
Get-Content .\.sync\sync.log -Tail 20     # what has it been doing?
Get-ScheduledTask BhuvanlabAutoSync        # is the watcher alive?
```

## Two deliberate safety behaviours

**1. The quiet period (45 s).** Regenerating the fourparam tables writes for
minutes. Committing mid-write would push a truncated CSV, so nothing is
committed until writes have actually stopped. Raise it for long batch runs:

```powershell
.\tools\install-autosync.ps1 -QuietSeconds 120
```

**2. The 250 MB push guard.** A single-tissue regeneration is ~42 MB and pushes
automatically. A *full* 100-table regeneration is ~1.9 GB raw and would saturate
the upstream link for a long time, so it is **committed locally** and the push
waits. Nothing is ever lost — release it deliberately:

```powershell
.\tools\sync-now.ps1 -Force
```

Changes to `fourparamsacrosstissues/data/` are held the same way, because those
are Git LFS pointers (see below) and unexpected changes there deserve a look.

## Repository facts worth knowing

**`data/` is Git LFS; this clone deliberately does not download it.** The 50
matrices are ~5.5 GB of LFS payload. This repo is configured with
`git lfs install --skip-smudge` and `lfs.fetchexclude=*`, so `data/*.csv.gz`
stay as small pointer files and no git operation will ever try to pull the
5.5 GB. To fetch one on purpose:

```powershell
git lfs pull --include="fourparamsacrosstissues/data/v11_log2_liver.csv.gz"
```

**`outputs/` is plain git, not LFS — on purpose.** These CSVs compress 2.46× in
git (1.94 GB → ~788 MB), while LFS stores blobs uncompressed and would consume
the full 1.94 GB of a quota that `data/` has already pushed to 5.5 GB. Plain git
is both smaller and cheaper here.

**This is a blobless partial clone.** It was attached with
`git fetch --filter=blob:none`, so `.git` is a few hundred KB instead of 813 MB —
the working files already matched the remote exactly, so downloading their blobs
would have been wasted. History is fetched on demand. If a command ever needs
old file contents it will fetch them automatically; to backfill everything:

```powershell
git fetch --refetch origin main
```

## Troubleshooting

**Nothing is syncing.** Check the task is running and read the log:

```powershell
Get-ScheduledTask BhuvanlabAutoSync | Select-Object State
Get-Content .\.sync\sync.log -Tail 30
```

**"Push HELD".** Expected for payloads over 250 MB or `data/` changes. Your work
is committed locally. Run `.\tools\sync-now.ps1 -Force`.

**"Diverged".** Someone (or another machine) pushed to `main` too. The script
refuses to guess — resolve it yourself:

```powershell
git pull --rebase
.\tools\sync-now.ps1
```

**Watch it work in the foreground** instead of as a hidden task:

```powershell
.\tools\install-autosync.ps1 -Uninstall
.\tools\sync-watch.ps1 -QuietSeconds 10
```
