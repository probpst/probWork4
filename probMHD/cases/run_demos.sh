#!/bin/bash
# turbulence-switch demonstration (3 runs, 1 rank each, sequential) then the 3-D torus (4 ranks)
set -e
cd "$(dirname "$0")"
B=$(pwd)/../probMHD
for t in 0 1 3; do
  d=KSTAR_turbulence_switch/turb$t
  rm -rf $d/Output_data $d/Volume_solutions_grids
  cp $B $d/
  (cd $d && mpirun --allow-run-as-root --oversubscribe -np 1 ./probMHD > log.txt 2>&1 && tail -2 log.txt | head -1)
done
(cd KSTAR_turbulence_switch && python3 compare.py)
d=KSTAR_3D_torus
rm -rf $d/Output_data $d/Volume_solutions_grids
cp $B $d/
(cd $d && mpirun --allow-run-as-root --oversubscribe -np 4 ./probMHD > log.txt 2>&1 && tail -2 log.txt | head -1 && python3 ../../tools/plot_monitor.py Output_data/KSTAR_3D_monitor.dat)
echo DEMOS_DONE
