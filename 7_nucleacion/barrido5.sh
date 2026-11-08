#!/bin/bash
# barrido5.sh — borde inferior del escalón Q = 1 (Claude, 2-11-2026)
# Pregunta: con λ = 24, ¿basta cualquier helicidad distinta de cero para atar un nudo de carga 1?
# PREDICCIÓN, escrita antes de ejecutarlo (con las cautelas de los barridos 2 y 3):
#   s = 0     → sin helicidad: Q = 0 (si apareciera algún nudo, su signo debería ser al azar)
#   s = 0,02 y 0,05 → no sé dónde está el borde; es justo lo que se mide
# Uso: bash barrido5.sh   (unos 30–40 minutos con 4 núcleos). Resultados en resultados5/.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
mkdir -p resultados5; cd resultados5; rm -f resumen.dat
COMUN="N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24"
for sg in 0 0.02 0.05; do for sd in 1 2 3; do
  ../nucleacion $COMUN s=$sg seed=$sd tag=s${sg}_sd$sd | grep FINAL; done; done
echo "== Resumen"; cat resumen.dat
