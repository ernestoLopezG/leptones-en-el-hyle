#!/bin/bash
# forma_capas_L155.sh — caja 15,5 (N = 200, mismo paso de malla) para ver si b/a sigue bajando como 1/L.
# Predicción anotada en resultados_forma.txt ANTES de ejecutarlo: b/a ≈ 1,23 en la capa n_z = 0.
# Uso: bash forma_capas_L155.sh   (del orden de 4–6 h con 4 núcleos; unos 200 MB de memoria)
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp ../7_nucleacion/nucleacion_pulsos.cpp -o nucleacion_pulsos
mkdir -p resultados_forma; cd resultados_forma
../nucleacion_pulsos N=200 L=15.5 T0=0 thold=0 tauq=0 trampa=0 lam=0 cada=200 cortes=1000000 init=hopf recupera=0 relax=3.0 tag=forma_L15.5 | grep FINAL
python3 ../forma_capas.py n_forma_L15.5.bin | tee forma_capas_L155.txt
