#!/bin/bash
# barrido_compacto.sh — comprobaciones de hopfion_compacto.cpp  (Claude, 4-11-2026)
# 1) caja mayor en x,y (L = 9): ¿el tamaño R_xy lo fija el nudo o la caja?
# 2) malla más fina (h ≈ 0,075): ¿el salto a «tubo enrollado» ocurre al mismo Lc? Si se mueve con h, es un artefacto de la malla.
# Cada etapa tarda unos minutos con 4 núcleos; el barrido completo, unas horas.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -fopenmp -o hopfion_compacto hopfion_compacto.cpp || exit 1
LCS="5.4 4.5 3.6 3.0 2.7 2.5 2.3 2.1 2.0 1.6 1.3 1.0"
./hopfion_compacto 101 60 9.0 1.0 10000 $LCS -t L9   | tee log_L9.txt
./hopfion_compacto 121 72 9.0 1.0 12000 $LCS -t L9f  | tee log_L9f.txt
python3 hopfq_compacto.py n_L9_*.bin n_L9f_*.bin | tee cargas.txt
python3 figuras_compacto.py compacto_L9.dat
python3 figuras_compacto.py compacto_L9f.dat
