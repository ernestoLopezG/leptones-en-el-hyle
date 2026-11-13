#!/bin/bash
# forma_capas.sh — forma de las capas del nudo relajado en función del tamaño de la caja (Claude, 13-11-2026)
# Mismo paso de malla h ≈ 0,078 en todas las cajas; relajación larga para acercarse al virial (E2/E4 → 1).
# En la nube (relax ≤ 1): b/a de las capas interiores ≈ 1,40–1,44 (caja 4,9), ≈ 1,35 (6,2), ≈ 1,30–1,31 (7,4):
# baja con la caja y no ha convergido. Este barrido dice hacia dónde va.
# Uso: bash forma_capas.sh   (del orden de 1,5–2 h con 4 núcleos)
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp ../7_nucleacion/nucleacion_pulsos.cpp -o nucleacion_pulsos
mkdir -p resultados_forma; cd resultados_forma
C="T0=0 thold=0 tauq=0 trampa=0 lam=0 cada=200 cortes=1000000 init=hopf recupera=0"
for NL in "80 6.2" "96 7.4" "128 9.9" "160 12.4"; do set -- $NL
  ../nucleacion_pulsos N=$1 L=$2 $C relax=3.0 tag=forma_L$2 | grep FINAL
done
python3 ../forma_capas.py n_forma_L6.2.bin n_forma_L7.4.bin n_forma_L9.9.bin n_forma_L12.4.bin | tee forma_capas.txt
