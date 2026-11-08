#!/bin/bash
# barrido.sh — barrido de nucleacion.cpp para el ordenador de casa (Claude, 2-11-2026)
# Uso:  bash barrido.sh            (todo; unas 2–3 horas con 4 núcleos)
#       bash barrido.sh rapido     (solo controles y 6 ejecuciones; unos 20 minutos)
# Resultados: resumen.dat (una línea por ejecución) y nucleacion_<tag>.png de cada una.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
echo "== Comprobaciones del programa"; ./nucleacion test
mkdir -p resultados; cd resultados; rm -f resumen.dat
COMUN="N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000"
corre() { ../nucleacion $COMUN "$@" | grep -E "^tag|FINAL"; python3 ../figuras_nucleacion.py "${@: -1}" >/dev/null 2>&1 || true; }

echo "== Control 1: un nudo ya formado sobrevive en esta malla y caja (debe quedar Q ≈ −1)"
../nucleacion $COMUN init=hopf T0=0 lam=0 thold=0 tauq=0 trampa=0 tag=ctl_hopf | grep FINAL
echo "== Control 2: fusión sin flujo (Kibble–Zurek puro); se espera Q = 0 casi siempre"
for sd in 1 2 3; do ../nucleacion $COMUN lam=0 seed=$sd tag=kz_s$sd | grep FINAL; done

if [ "$1" == "rapido" ]; then
  for sg in 1 -1; do for sd in 1 2 3; do corre lam=20 s=$sg seed=$sd tag=l20_s${sg}_sd$sd; done; done
else
  # Pregunta 1: ¿el signo de la helicidad elige el signo de Q?   Pregunta 2: ¿qué λ da |Q| = 1 y no 2–3?
  # Pregunta 3: ¿qué fracción sobrevive al retirar el flujo?
  for l in 12 16 20 24; do for sg in 1 -1; do for sd in 1 2 3; do
      corre lam=$l s=$sg seed=$sd tag=l${l}_s${sg}_sd$sd; done; done; done
  # Anillo sin remolino (helicidad 0): ¿se forma algo?
  for sd in 1 2 3; do corre lam=20 s=0 seed=$sd tag=l20_s0_sd$sd; done
  # Flujo permanente (la pulsación sigue ahí): ¿el nudo sostenido es estable?
  for sd in 1 2 3; do corre lam=16 s=1 retira=0 seed=$sd tag=l16_perm_sd$sd; done
  # Enfriamiento más rápido (congela antes)
  for sd in 1 2 3; do corre lam=16 s=1 tauq=0.05 seed=$sd tag=l16_rapido_sd$sd; done
fi
python3 ../figuras_nucleacion.py resumen >/dev/null 2>&1 || true
echo "== Resumen"; cat resumen.dat
