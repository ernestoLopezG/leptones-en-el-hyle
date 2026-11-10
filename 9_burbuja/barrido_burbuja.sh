#!/bin/bash
# barrido_burbuja.sh — ¿sobrevive el nudo en una burbuja ordenada dentro de un fondo desordenado?  (Claude, 10-11-2026)
# Malla h ≈ 0,09 (la de la tabla de convergencia; más gruesa, el nudo se deshace en la malla).
# Tiempo de relajación τ = 0,02·h²·iteraciones ≈ 6,5 en todas las corridas (40 000 iteraciones).
# Unas 4,5–5,5 horas con 4 núcleos (i5-4590); la más larga, Rb = 8 (casi 2 h), necesita unos 2,5 GB de memoria.
# Los ficheros burbuja_*.dat guardan la energía cada 200 iteraciones: aunque se corte, lo hecho sirve.
# Uso: bash barrido_burbuja.sh        Resultados: resultados_burbuja/resumen_burbuja.dat y q_burbuja.txt
set -e
g++ -O3 -march=native -fopenmp hopfion_burbuja.cpp -o hopfion_burbuja
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
mkdir -p resultados_burbuja && cd resultados_burbuja
IT=40000; W=0.4
corre () {   # Rb fondo w etiqueta
  local Rb=$1 fondo=$2 w=$3 et=$4
  local L=$(python3 -c "print(round(2*($Rb + 5*$w + 0.0), 3))")
  local N=$(python3 -c "print(int(round($L/0.09)) + 1)")
  echo "== $et: Rb=$Rb w=$w fondo=$fondo L=$L N=$N  ($(date +%H:%M))"
  ../hopfion_burbuja $N $L 1.0 $Rb $w $fondo $IT $et > log_$et.txt
  tail -1 log_$et.txt
  python3 ../hopfq_burbuja.py n_$et.bin | tee -a q_burbuja.txt
}
# 1) Fondo uniforme fuera de la burbuja (sin orden pero con el mismo valor que el nudo en el infinito)
for Rb in 2 3 4 6 8; do corre $Rb uniforme $W u_Rb$Rb; done
# 2) Fondo aleatorio (direcciones al azar fuera de la burbuja): el caso del hyle
for Rb in 3 4 6; do corre $Rb aleatorio $W a_Rb$Rb; done
# 3) Sensibilidad a la anchura del borde
corre 3 uniforme 0.8 u_Rb3_w08
echo "Terminado $(date +%H:%M). Resumen:"; cat resumen_burbuja.dat; cat q_burbuja.txt
