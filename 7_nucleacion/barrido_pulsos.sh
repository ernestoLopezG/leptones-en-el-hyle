#!/bin/bash
# barrido_pulsos.sh — ¿sobrevive el nudo a una fusión PULSADA?  (Claude, 7-11-2026)
#
# Hipótesis de Ernesto (7-11): el modo 3/2, más intenso, produce nucleaciones permanentes que desatan su propio
# nudo; el ½ solo produce nucleaciones breves. Bloque A5: un núcleo fundido permanente sobre el anillo del nudo
# lo desata incluso con κ₄. Aquí el núcleo se enciende una fracción frac de cada periodo.
#
# CALIBRACIÓN hecha en la nube (malla 64, caja 4,9, Tp = 1,5, núcleo toro sobre el anillo, theta = 0,05):
#   con fusión permanente el nudo se desata en t ≈ 0,02–0,05;
#   periodo 0,004: frac 0,25 sobrevive, 0,5 se desata     periodo 0,01: frac 0,25 sobrevive, 0,5 se desata
#   periodo 0,04:  frac 0,10 sobrevive, 0,25 se desata    periodo 0,1:  frac 0,10 se desata (encendido 0,01)
#   → umbral doble: fracción de tiempo fundido entre 0,25 y 0,5, y tiempo encendido por pulso ≲ 0,005.
# Este barrido afina el umbral en la malla de producción y comprueba de qué depende (ruido, semilla).
#
# Uso: bash barrido_pulsos.sh   (unos 30 minutos con 4 núcleos). Resultados en resultados_pulsos/.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion_pulsos.cpp -o nucleacion_pulsos
mkdir -p resultados_pulsos; cd resultados_pulsos; rm -f resumen.dat
P=../nucleacion_pulsos
C="N=80 L=6.2 T0=0 thold=0 tauq=0 trampa=0 lam=0 cada=50 cortes=100000"
NUC="Tp=1.5 forma=toro Rt=1.0 Rp=0.35 relax=0.3 recupera=0.1"
c() { $P $C "$@" | grep FINAL; }
[ -f n_base.bin ] || c init=hopf relax=0.5 recupera=0 tag=base
echo "== 1) umbral en la fracción de tiempo fundido (pulsos rápidos, periodo 0,004)"
for f in 0.2 0.3 0.35 0.4 0.45; do c init=n_base.bin $NUC periodo=0.004 frac=$f tag=P0.004_f$f; done
echo "== 2) umbral en la duración de cada pulso (fracción fija 0,2)"
for p in 0.01 0.02 0.03 0.05; do c init=n_base.bin $NUC periodo=$p frac=0.2 tag=P${p}_f0.2; done
echo "== 3) ¿depende del ruido térmico? (periodo 0,004, fracción 0,3)"
for th in 0 0.025 0.1; do c init=n_base.bin $NUC periodo=0.004 frac=0.3 theta=$th tag=th${th}_P0.004_f0.3; done
echo "== 4) control: fusión permanente SIN ruido (¿basta perder la rigidez, o hace falta el ruido?)"
c init=n_base.bin $NUC theta=0 tag=perm_th0
echo "== 5) semillas (periodo 0,004, fracción 0,3)"
for sd in 2 3; do c init=n_base.bin $NUC periodo=0.004 frac=0.3 seed=$sd tag=sd${sd}_P0.004_f0.3; done
echo "== Resumen"; cat resumen.dat
