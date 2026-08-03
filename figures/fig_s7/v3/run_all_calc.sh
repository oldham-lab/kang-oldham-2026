#!/usr/bin/env bash
# Recompute per-celltype mean/SE for all seven Fig. S7 series, restricted to the
# common gene universe (v3/common_genes.txt). Runs the four dataset scripts in
# parallel and waits for all to finish. Per-script logs in v3/logs/.
set -u
cd "$(dirname "$0")"
mkdir -p logs
pids=()
for s in jorstad seaad MIT_AD morabito; do
  echo "launching $s ..."
  python3 "calc_mean_se_per_celltype_${s}.py" > "logs/${s}.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do
  wait "$p" || fail=1
done
echo "=== all calc jobs done (fail=$fail) ==="
tail -n 3 logs/*.log
exit $fail
