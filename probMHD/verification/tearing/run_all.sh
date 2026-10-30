#!/bin/bash
# runs the cases two at a time (2 MPI ranks each)
cd "$(dirname "$0")"
python3 make_cases.py
run() { cp ../../probMHD "$1"/; (cd "$1" && mpirun --allow-run-as-root --oversubscribe -np 2 ./probMHD > log.txt 2>&1; tail -2 log.txt | head -1); }
run ideal_28 & run ideal_56 & wait
run S100_28 & run S100_56 & wait
run S300_56 & run S300_56_eno3 & wait
python3 plot_tearing.py
echo TEARING_DONE
