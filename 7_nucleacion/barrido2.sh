#!/bin/bash
# barrido2.sh — segundo barrido de nucleacion.cpp (Claude, 2-11-2026)
# Pregunta 1: ¿con menos remolino (|s| < 1) el nudo atado es de carga 1 en lugar de 3?
# Pregunta 2: la asimetría entre s = +1 y s = −1 del primer barrido, ¿es azar? (por simetría deberían ser
#             equivalentes con Q → −Q; se repite s = ±1 con más semillas)
# Uso: bash barrido2.sh     (unas 2–3 horas con 4 núcleos). Resultados en resultados2/.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
mkdir -p resultados2; cd resultados2; rm -f resumen.dat
COMUN="N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24"
corre() { ../nucleacion $COMUN "$@" | grep -E "FINAL"; python3 ../figuras_nucleacion.py "${@: -1}" >/dev/null 2>&1 || true; }
for sg in 0.3 -0.3 0.5 -0.5 0.7 -0.7; do for sd in 1 2 3 4; do corre s=$sg seed=$sd tag=s${sg}_sd$sd; done; done
for sg in 1 -1; do for sd in 4 5 6; do corre s=$sg seed=$sd tag=s${sg}_sd$sd; done; done
python3 ../figuras_nucleacion.py resumen >/dev/null 2>&1 || true
echo "== Resumen"; cat resumen.dat
