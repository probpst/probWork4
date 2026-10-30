#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 make_cases.py
for d in S300 S1000; do
  cp ../../probMHD "$d"/
  (cd "$d" && mpirun --allow-run-as-root --oversubscribe -np 2 ./probMHD > log.txt 2>&1 && tail -2 log.txt | head -1)
done
python3 plot_tearing.py
