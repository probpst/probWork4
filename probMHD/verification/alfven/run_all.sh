#!/bin/bash
# run every case folder with the solver copied from the package root
set -e
cd "$(dirname "$0")"
python3 make_cases.py
for d in N*_s*; do
  cp ../../probMHD "$d"/
  (cd "$d" && mpirun --allow-run-as-root --oversubscribe -np 1 ./probMHD > log.txt 2>&1 && tail -2 log.txt | head -1)
done
python3 plot_alfven.py
