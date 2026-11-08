#!/bin/bash
# convergencia.sh — tabla de convergencia de la sección 4 (malla y caja)  (Claude, 7-11-2026)
# Uso: bash convergencia.sh   (unas 2 horas con 2 núcleos; la malla 121 necesita continuar, ver abajo)
set -e
g++ -O3 -march=native -fopenmp hopfion.cpp -o hopfion
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
for cfg in "61 5.4" "81 7.2" "101 9.0" "121 9.0"; do set -- $cfg
  mkdir -p N$1_L$2 && cd N$1_L$2
  ../hopfion $1 $2 1.0 20000 > log.txt 2>&1
  python3 ../hopfq.py n_final_N$1.bin > q.txt; python3 ../medir_tamano.py n_final_N$1.bin > tam.txt; cd ..
  echo "N=$1 L=$2: $(tail -1 N$1_L$2/hopf_N$1.dat)"
done
# Con dt ∝ h², la malla 121 relaja más despacio: se continúa 16 000 iteraciones más.
cd N121_L9.0 && ../hopfion 121 9.0 1.0 16000 n_final_N121.bin > log2.txt 2>&1 && tail -1 hopf_N121.dat && cd ..
# Prueba de punto de silla (sección 4, «Estabilidad»): deformar el nudo de la malla 61 y relajarlo otra vez.
cd N61_L5.4 && python3 ../perturba.py n_final_N61.bin && mkdir -p silla && mv perturbado.bin silla/ && cd silla \
  && ../../hopfion 61 5.4 1.0 20000 perturbado.bin > log.txt 2>&1 && tail -1 hopf_N61.dat && python3 ../../perturba.py --asimetria n_final_N61.bin
