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
#      bash barrido_ablandamiento.sh solo4    (solo la formación, bloque 4: unas 1–1,5 h con 4 núcleos)
#      bash barrido_ablandamiento.sh solo5    (¿el Q = 2 sin fusión es artefacto?: unas 1,5–2 h)
#      bash barrido_ablandamiento.sh solo6    (control: ¿la fusión sigue atando en N = 112?: unos 30–40 min)
# Resultados: resultados_ablandamiento/resumen.dat (una línea por ejecución) y salida en pantalla.
set -e
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
g++ -O3 -march=native -fopenmp nucleacion_pulsos.cpp -o nucleacion_pulsos
g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
mkdir -p resultados_ablandamiento; cd resultados_ablandamiento; [ "$1" = "solo4" ] || [ "$1" = "solo5" ] || [ "$1" = "solo6" ] || rm -f resumen.dat
P=../nucleacion_pulsos; F=../nucleacion
NUC="forma=toro Rt=1.0 Rp=0.35 relax=0.3 recupera=0.1"
c() { $P "$@" | grep FINAL; }

if [ "$1" = "solo6" ]; then
echo "########## 6) control: fusión (T0 = 2) en N = 112, caja 7,4, s = 0,4 ##########"
# En N = 112 el ablandamiento T0 = 0,9 no ató nada. Falta saber si la fusión sí ata en esa malla;
# si no atara, el problema sería de la malla fina en general, no del ablandamiento.
$F N=112 L=7.4 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24 s=0.4 T0=2.0 seed=1 tag=f6_N112_T02.0_s0.4 | grep FINAL
$F N=112 L=7.4 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24 s=0.4 T0=0.9 seed=2 tag=f6_N112_T00.9_s0.4_sd2 | grep FINAL
echo "== Resumen"; cat resumen.dat; exit 0
fi

if [ "$1" = "solo5" ]; then
echo "########## 5) ¿el Q = 2 sin fusión (T0 = 0,9, s = 0,4) es artefacto de malla? ##########"
# En N = 80 el ablandamiento profundo con flujo muy retorcido ató un nudo de carga 2 idéntico al de la fusión
# (E2+E4 = 924,8). Sin fusión la topología lo prohíbe en el continuo. Pruebas: malla más fina con el mismo h
# relativo (N = 96, caja 7,4, y N = 112, caja 7,4) y ablandamientos menos profundos.
for t0 in 2.0 0.9; do $F N=96 L=7.4 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24 s=0.4 T0=$t0 seed=1 tag=f5_N96_T0${t0}_s0.4 | grep FINAL; done
$F N=112 L=7.4 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24 s=0.4 T0=0.9 seed=1 tag=f5_N112_T00.9_s0.4 | grep FINAL
for t0 in 0.7 0.5 0.3; do $F N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24 s=0.4 T0=$t0 seed=1 tag=f5_N80_T0${t0}_s0.4 | grep FINAL; done
echo "== Resumen"; cat resumen.dat; exit 0
fi
if [ "$1" != "solo4" ]; then
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

fi
[ "$1" = "rapido" ] && { echo "== Resumen"; cat resumen.dat; exit 0; }

echo "########## 4) formación desde el vacío: fusión (T0=2) frente a ablandamiento (T0=0,9) ##########"
# Configuración de producción de barrido3.sh (allí la fusión con s = 0,2–0,4 ató nudos de carga 1): caja 6,2,
# anillo de flujo Rc = 1,8. La versión anterior de este bloque usaba caja 4,9, demasiado pequeña para ese anillo,
# y con fusión tampoco se formaba nada (sin control positivo el bloque no informaba); en N = 96 con h = 0,05 el
# esquema explícito se volvía inestable (E2+E4 ~ 3·10^5).
COMUN="N=80 L=6.2 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24"
for t0 in 2.0 0.9; do for sg in 0.2 0.4; do for sd in 1 2; do
  $F $COMUN s=$sg T0=$t0 seed=$sd tag=form_T0${t0}_s${sg}_sd$sd | grep FINAL
done; done; done
# comprobación de malla con el mismo paso h que en producción (N=96, caja 7,4)
for t0 in 2.0 0.9; do $F N=96 L=7.4 Rc=1.8 relax=0.5 cada=100 cortes=2000 lam=24 s=0.2 T0=$t0 seed=1 tag=form96_T0${t0}_s0.2 | grep FINAL; done

echo "== Resumen"; cat resumen.dat
