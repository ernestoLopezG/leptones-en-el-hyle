#!/bin/bash
# barrido_persistente.sh — nucleación persistente y forzado ½ frente a 3/2  (Claude, 7-11-2026)
#
# Tres preguntas (objeciones M3 y M5 del árbitro; respuestas_arbitro.md):
#   A) ¿Puede un núcleo fundido permanente sustituir al término κ₄ (E4) para que el nudo no se encoja?
#   B) Desde el vacío y con la MISMA energía de forzado, ¿ata un nudo el forzado «½» y no el «3/2»?
#   C) Con un nudo Q = 1 ya formado, ¿puede el forzado «3/2» convertirlo en Q = 3? ¿Lo impide o lo facilita
#      la nucleación persistente?
#
# Correspondencia usada [H]: el director es bilineal en el espinor (k_director = 2·k_espinor), así que un modo
# de orden ½ corresponde a un nudo Q = 1 y uno de orden 3/2 a Q = 3. En los barridos 2–5 el nudo atado sigue
# la torsión del flujo, |Q| ≈ redondeo de 2T, y para w = 0,6 es 2T = 3,189·s (torsion_flujo.py). Por tanto:
#   forzado ½:   s = 0,3136 (2T = 1),  λ = 24
#   forzado 3/2: s = 0,9407 (2T = 3),  λ = 18,32   ← misma energía de forzado: λ²(1 + s²) igual en los dos
#
# PREDICCIONES, escritas ANTES de ejecutarlo (Claude):
#   A1 (k4=1, sin núcleo): el nudo sobrevive, Q ≈ −1 (control).
#   A2 (k4=0, sin núcleo): Derrick: se encoge hasta la malla y Q → 0.
#   A3/A4 (k4=0, núcleo fundido toro/bola): creo que FALLA. En la zona fundida el director no tiene rigidez,
#      así que es un agujero por el que el nudo se puede desatar. Una prueba de humo (malla 48, t = 0,05) ya
#      mostró Q: −0,93 → −0,46. Si falla, «núcleo fundido en lugar de κ₄» queda refutado en este juguete.
#   A5 (k4=1, núcleo toro): puede degradar incluso un nudo bueno, por la misma razón.
#   B: con la misma energía, ½ → |Q| = 1 y 3/2 → |Q| = 3 (como en los barridos 2–5). Si es así, el mecanismo
#      (2) de Ernesto, «formar e⁻e⁺e⁻ es más fácil que el estado excitado», NO recibe apoyo de este juguete.
#   C: sin núcleo (Tp=0), Q se queda en −1: en un medio ordenado Q no puede cambiar sin fundirse (la barrera es
#      el orden, no la nucleación). Con núcleo persistente, la puerta está abierta y Q PUEDE saltar a −3.
#      Es decir: espero lo contrario de la hipótesis (1) en su forma literal. Si el núcleo persistente impidiera
#      el salto, la hipótesis (1) tendría apoyo y mi predicción sería errónea.
#
# Uso: bash barrido_persistente.sh        (unas 2,5 horas con 4 núcleos). Resultados en resultados_persistente/.
#      bash barrido_persistente.sh A      (solo el bloque A, unos 30 min)
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion_persistente.cpp -o nucleacion_persistente
./nucleacion_persistente test
mkdir -p resultados_persistente; cd resultados_persistente
P=../nucleacion_persistente
MALLA="N=80 L=6.2 cada=100 cortes=4000"
QUIETO="T0=0 thold=0 tauq=0 trampa=0"          # sin fusión inicial: empieza directamente en la fase 3
c() { $P $MALLA "$@" | grep FINAL; }
BLOQUE=${1:-ABC}

# 0) nudo Q = 1 relajado, punto de partida de A y C
[ -f n_base.bin ] || c init=hopf $QUIETO lam=0 relax=1 tag=base

if [[ $BLOQUE == *A* ]]; then echo "== A: núcleo fundido en lugar de κ₄"
  c init=n_base.bin $QUIETO lam=0 relax=1 k4=1                             tag=A1_k4_sin
  c init=n_base.bin $QUIETO lam=0 relax=1 k4=0                             tag=A2_sinK4_sin
  c init=n_base.bin $QUIETO lam=0 relax=1 k4=0 Tp=1.5 forma=toro Rt=1.0 Rp=0.35 tag=A3_sinK4_toro
  c init=n_base.bin $QUIETO lam=0 relax=1 k4=0 Tp=1.5 forma=bola Rp=0.8    tag=A4_sinK4_bola
  c init=n_base.bin $QUIETO lam=0 relax=1 k4=1 Tp=1.5 forma=toro Rt=1.0 Rp=0.35 tag=A5_k4_toro
fi

if [[ $BLOQUE == *B* ]]; then echo "== B: desde el vacío, forzado ½ frente a 3/2 con la misma energía (flujo permanente)"
  VAC="Rc=1.8 retira=0 relax=0.5"
  for sd in 1 2 3; do
    c $VAC s=0.3136 lam=24    seed=$sd                          tag=B_medio_sd$sd
    c $VAC s=0.9407 lam=18.32 seed=$sd                          tag=B_3medios_sd$sd
    c $VAC s=0.3136 lam=24    seed=$sd Tp=1.5 forma=bola Rp=0.6 tag=B_medio_per_sd$sd
    c $VAC s=0.9407 lam=18.32 seed=$sd Tp=1.5 forma=bola Rp=0.6 tag=B_3medios_per_sd$sd
  done
fi

if [[ $BLOQUE == *C* ]]; then echo "== C: nudo Q = 1 sometido a forzado 3/2 (flujo permanente)"
  c init=n_base.bin $QUIETO retira=0 relax=1.5 s=0.3136 lam=24                       tag=C_medio_control
  c init=n_base.bin $QUIETO retira=0 relax=1.5 s=0.9407 lam=18.32                    tag=C_3medios_sin
  for sd in 1 2; do
    c init=n_base.bin $QUIETO retira=0 relax=1.5 s=0.9407 lam=18.32 seed=$sd Tp=1.5 forma=toro Rt=1.0 Rp=0.35 tag=C_3medios_toro_sd$sd
    c init=n_base.bin $QUIETO retira=0 relax=1.5 s=0.9407 lam=18.32 seed=$sd Tp=1.5 forma=bola Rp=0.6 tag=C_3medios_bola_sd$sd
  done
fi
echo "== Resumen"; cat resumen.dat
for f in log_*.dat; do python3 ../figuras_persistente.py "$f" >/dev/null 2>&1 || true; done
