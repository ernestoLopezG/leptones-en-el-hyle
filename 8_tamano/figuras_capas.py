#!/usr/bin/env python3
# figuras_capas.py — resumen del barrido de hopfion_capas.cpp  (Claude, 4-11-2026)
# Lee resumen_capas.dat. Para cada corrida con capas, divide R_xy entre el de su control (mismo Delta, B = 0)
# y lo dibuja frente a √(λΔ). Si la elasticidad de las capas fijara el tamaño según √(λΔ), los puntos
# seguirían una recta por el origen; si no, R_xy/R_control queda independiente de λ y Δ.
import numpy as np, collections
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
filas = [l.split() for l in open('resumen_capas.dat') if l.strip()]
d = collections.OrderedDict()
for f in filas: d[f[0]] = dict(D=float(f[1]), B=float(f[2]), lam=float(f[3]), borde=f[5], Efs=float(f[6]),
                              Ec=float(f[7]), R=float(f[8]), sz=float(f[9]), sl=float(f[10]), it=int(f[11]))
ctrl = {v['D']: v['R'] for v in d.values() if v['B'] == 0}
print(f"{'corrida':22s} {'Δ':>4s} {'B':>4s} {'λ':>4s} {'borde':>5s} {'R_xy':>6s} {'R/Rctrl':>7s} {'σ_z':>6s} {'√(λΔ)':>6s} {'E_FS':>7s} {'E_capas':>7s} {'it':>6s}")
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for k, v in d.items():
    rr = v['R']/ctrl[v['D']] if v['D'] in ctrl else np.nan
    print(f"{k:22s} {v['D']:4.1f} {v['B']:4g} {v['lam']:4g} {v['borde']:>5s} {v['R']:6.3f} {rr:7.3f} {v['sz']:6.3f} {v['sl']:6.3f} {v['Efs']:7.4f} {v['Ec']:7.4f} {v['it']:6d}")
    if v['B'] > 0:
        m = 'o' if v['borde'] == 'fijo' else 's'
        ax[0].scatter(v['sl'], rr, marker=m, s=20 + 8*v['B'], label=k)
        ax[1].scatter(v['B'], rr, marker=m, label=k)
ax[0].set_xlabel('√(λΔ)'); ax[0].set_ylabel('R_xy / R_xy(control, mismo Δ)'); ax[0].grid(alpha=.3)
ax[0].set_title('¿Sigue el tamaño extendido a √(λΔ)?'); ax[0].legend(fontsize=7)
ax[1].set_xscale('log'); ax[1].set_xlabel('B (módulo de compresión)'); ax[1].set_ylabel('R_xy / R_control'); ax[1].grid(alpha=.3)
ax[1].set_title('Efecto de la rigidez de las capas')
plt.tight_layout(); plt.savefig('capas.png', dpi=120); print('-> capas.png')
