#!/bin/bash
# barrido3.sh — bordes del escalón Q = 1 (Claude, 2-11-2026)
# Predicción hecha ANTES de ejecutarlo (con torsion_flujo.py, si Q ≈ redondeo de 2T):
#   s = 0,10 → 2T = 0,32 → Q = 0      s = 0,15 → 2T = 0,48 → frontera 0/1
#   s = 0,20 → 2T = 0,64 → Q = 1      s = 0,40 → 2T = 1,27 → Q = 1 (quizá alguna vez 2)
# Uso: bash barrido3.sh   (unas 1,5 horas con 4 núcleos). Resultados en resultados3/.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
mkdir -p resultados3; cd resultados3; rm -f resumen.dat
COMUN="N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24"
for sg in 0.1 0.15 0.2 0.4; do for sd in 1 2 3; do
  ../nucleacion $COMUN s=$sg seed=$sd tag=s${sg}_sd$sd | grep FINAL; done; done
for sd in 1 2; do ../nucleacion $COMUN s=-0.2 seed=$sd tag=s-0.2_sd$sd | grep FINAL; done
echo "== Resumen"; cat resumen.dat
