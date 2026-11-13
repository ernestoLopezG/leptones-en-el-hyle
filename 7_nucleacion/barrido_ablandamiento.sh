#!/bin/bash
# barrido_ablandamiento.sh — ablandamiento (T < 1) frente a fusión (T > 1) en mallas finas  (Claude, 13-11-2026)
#
# Pregunta: las zonas de mayor índice son ablandamientos del orden, no fusiones (S² = 1 − T > 0). ¿Qué le hacen al
# nudo? Las pruebas rápidas en malla gruesa (N=64, resultados_ablandamiento.txt) dieron:
#   - ablandamiento a pulsos sobre el anillo: el nudo sobrevive intacto (la fusión con el mismo pulso lo desata);
#   - ablandamiento permanente: lo desata incluso sin ruido (T ≥ 0,6), probablemente porque el término de cuarto
#     orden va pesado con S⁶ y el de segundo con S², y el nudo pierde la protección frente a Derrick;
#   - sin fusión no se forma ningún nudo desde el vacío (el Q = −2,9 de N=64 desaparece en N=80: artefacto).
# Este barrido lo comprueba en mallas más finas y con más semillas.
#
# Uso: bash barrido_ablandamiento.sh          (todo: del orden de 1,5–2 h con 4 núcleos)
#      bash barrido_ablandamiento.sh rapido   (solo los bloques 1–3 en N=80: unos 20 min)
# Resultados: resultados_ablandamiento/resumen.dat (una línea por ejecución) y salida en pantalla.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion_pulsos.cpp -o nucleacion_pulsos
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
mkdir -p resultados_ablandamiento; cd resultados_ablandamiento; rm -f resumen.dat
P=../nucleacion_pulsos; F=../nucleacion
NUC="forma=toro Rt=1.0 Rp=0.35 relax=0.3 recupera=0.1"
c() { $P "$@" | grep FINAL; }

for N in 80 96; do
  [ "$1" = "rapido" ] && [ $N = 96 ] && break
  C="N=$N L=6.2 T0=0 thold=0 tauq=0 trampa=0 lam=0 cada=100 cortes=1000000"
  echo "########## malla N=$N, caja 6,2 ##########"
  [ -f n_base_N$N.bin ] || $P $C init=hopf relax=0.5 recupera=0 tag=base_N$N | grep FINAL
  B=n_base_N$N.bin

  echo "== 1) ablandamiento PERMANENTE sobre el anillo, sin ruido (control: fusión Tp=1,5)"
  for tp in 1.5 0.9 0.7 0.5 0.3; do c $C init=$B $NUC Tp=$tp theta=0 tag=N${N}_perm_Tp$tp; done

  echo "== 2) ablandamiento A PULSOS (periodo 0,004 y 0,04; mitad del tiempo encendido), sin ruido"
  for tp in 1.5 0.9 0.6; do for p in 0.004 0.04; do
    c $C init=$B $NUC Tp=$tp theta=0 periodo=$p frac=0.5 tag=N${N}_pul_Tp${tp}_P$p; done; done

  echo "== 3) ¿el ablandamiento permanente débil (Tp=0,3) acaba desatándolo? (relajación larga)"
  c $C init=$B $NUC Tp=0.3 theta=0 relax=1.5 tag=N${N}_perm_Tp0.3_largo
done

[ "$1" = "rapido" ] && { echo "== Resumen"; cat resumen.dat; exit 0; }

echo "########## 4) formación desde el vacío: fusión (T0=2) frente a ablandamiento (T0=0,9), flujo intenso ##########"
for N in 80 96; do for t0 in 2.0 0.9; do for sd in 1 2; do
  $F N=$N L=4.9 Rc=1.8 lam=24 s=1 T0=$t0 seed=$sd relax=0.3 cada=200 cortes=1000000 tag=form_N${N}_T0${t0}_sd$sd | grep FINAL
done; done; done

echo "== Resumen"; cat resumen.dat
