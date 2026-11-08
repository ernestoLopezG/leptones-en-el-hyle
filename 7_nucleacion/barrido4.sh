#!/bin/bash
# barrido4.sh — ¿el factor 2 de Q ≈ 2T es real o es la geometría del tubo del flujo? (Claude, 2-11-2026)
# Se cambia el grosor del tubo (w = 0,4 y 0,8; en los barridos anteriores w = 0,6) y se elige s para que la
# torsión del flujo sea la misma (calculada con torsion_flujo.py):
#     w = 0,4: T = 2,729·s  →  T = 0,5 con s = 0,183;  T = 1,0 con s = 0,366
#     w = 0,8: T = 1,011·s  →  T = 0,5 con s = 0,495;  T = 1,0 con s = 0,989
# La intensidad del flujo sobre el nudo crece con el área del tubo (∝ w²). Se prueban dos opciones:
#     λ = 24 (como antes) y λ reescalado para el mismo empuje total: λ·w² = 24·0,36  →  λ = 54 (w = 0,4), 13,5 (w = 0,8)
# PREDICCIÓN, escrita antes de ejecutarlo: si el factor 2 es real, T = 0,5 → Q = 1 y T = 1,0 → Q = 2
#   en los dos grosores. Si Q cambia con w a torsión fija, el factor 2 venía de la geometría.
# Uso: bash barrido4.sh   (algo más de 1 hora con 4 núcleos). Resultados en resultados4/.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
mkdir -p resultados4; cd resultados4; rm -f resumen.dat
COMUN="N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000"
c() { ../nucleacion $COMUN "$@" | grep FINAL; }
for sd in 1 2; do
  c w=0.4 s=0.183 lam=24   seed=$sd tag=w0.4_T0.5_l24_sd$sd
  c w=0.4 s=0.183 lam=54   seed=$sd tag=w0.4_T0.5_l54_sd$sd
  c w=0.4 s=0.366 lam=54   seed=$sd tag=w0.4_T1.0_l54_sd$sd
  c w=0.8 s=0.495 lam=24   seed=$sd tag=w0.8_T0.5_l24_sd$sd
  c w=0.8 s=0.495 lam=13.5 seed=$sd tag=w0.8_T0.5_l13.5_sd$sd
  c w=0.8 s=0.989 lam=13.5 seed=$sd tag=w0.8_T1.0_l13.5_sd$sd
done
echo "== Resumen"; cat resumen.dat
