#!/bin/bash
# barrido_capas.sh — ¿fija la elasticidad de las capas el tamaño del nudo?  (Claude, 4-11-2026)
# Malla h = 0,09 en las tres direcciones (con h ≥ 0,1 el nudo se deshace en la malla).
# Cada corrida, hasta 30 000 iteraciones: de 10 a 20 min con 4 núcleos. Total, unas 3-4 horas.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -fopenmp -o hopfion_capas hopfion_capas.cpp -lfftw3 -lm || exit 1
./hopfion_capas test || exit 1
IT=30000
c() { ./hopfion_capas 101 $1 9 $2 1.0 $3 $4 1 $5 $IT $6; }
#   Nz  Delta  B   lambda borde etiqueta
# A) dependencia con lambda y con B (Delta = 4,5, capas fijas en las dos interfaces)
c 51 4.5  0   0.3  fijo  D4.5_control
c 51 4.5  4   0.1  fijo  D4.5_B4_l0.1
c 51 4.5  4   0.3  fijo  D4.5_B4_l0.3
c 51 4.5  4   1.0  fijo  D4.5_B4_l1.0
c 51 4.5  1   0.3  fijo  D4.5_B1_l0.3
c 51 4.5  16  0.3  fijo  D4.5_B16_l0.3
# B) dependencia con el espesor Delta
c 41 3.6  0   0.3  fijo  D3.6_control
c 41 3.6  4   0.3  fijo  D3.6_B4_l0.3
c 61 5.4  0   0.3  fijo  D5.4_control
c 61 5.4  4   0.3  fijo  D5.4_B4_l0.3
# C) capas libres en la interfaz exterior (sin la frustración del modo uniforme)
c 41 3.6  4   0.3  libre D3.6_B4_l0.3_libre
c 51 4.5  4   0.3  libre D4.5_B4_l0.3_libre
c 61 5.4  4   0.3  libre D5.4_B4_l0.3_libre
python3 figuras_capas.py
