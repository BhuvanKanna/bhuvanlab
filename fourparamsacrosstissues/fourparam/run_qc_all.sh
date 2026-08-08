#!/usr/bin/env bash
# Build QC tables for the top-N tissues by donor count, raw and excluded.
#
# Safe to re-run at any time: compute_qc.py skips a table whose final CSV
# exists, and resumes from `<out>.partial` otherwise, so a killed run loses at
# most one flush interval. An output_lock stops two copies of this script from
# racing on the same table.
#
#   ./run_qc_all.sh                 # all tissues below, both kinds
#   ./run_qc_all.sh muscle_skeletal # just one
#
# Progress goes to stdout; each table also prints its own dist_class census.
set -u

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

JOBS="${JOBS:-10}"          # 12 cores here; leave a couple for everything else
FLUSH="${FLUSH:-1000}"

TISSUES=(
  muscle_skeletal                    # 818 donors
  whole_blood                        # 803
  skin_sun_exposed_lower_leg         # 754
  adipose_subcutaneous               # 714
  artery_tibial                      # 691
  thyroid                            # 684
  nerve_tibial                       # 670
  cells_cultured_fibroblasts         # 652
  skin_not_sun_exposed_suprapubic    # 651
  esophagus_mucosa                   # 614
)
[ $# -gt 0 ] && TISSUES=("$@")

for t in "${TISSUES[@]}"; do
  m="../data/v11_log2_${t}.csv.gz"
  if [ ! -f "$m" ]; then
    echo "SKIP $t: matrix missing"; continue
  fi
  if [ "$(stat -c%s "$m")" -lt 1000 ]; then
    echo "SKIP $t: matrix is an LFS pointer, not the data"; continue
  fi
  for kind in raw excluded; do
    if [ "$kind" = raw ]; then THR=(); else THR=(--threshold -1); fi
    echo "=============== $t / $kind  ($(date +%H:%M:%S)) ==============="
    python -u compute_qc.py --input "$m" --id-col Name --name-col Description \
        --jobs "$JOBS" --flush-every "$FLUSH" "${THR[@]}" || echo "FAILED: $t/$kind"
  done
done

echo "=============== done ($(date +%H:%M:%S)) ==============="
ls -1 ../qc/*.csv 2>/dev/null | wc -l | xargs echo "qc tables present:"
