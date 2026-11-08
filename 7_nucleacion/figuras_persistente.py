#!/usr/bin/env python3
# figuras_persistente.py — Q, Q_S, tamaño y energía frente al tiempo para un log de nucleacion_persistente
# Uso: python3 figuras_persistente.py log_<tag>.dat   → persistente_<tag>.png
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
f = sys.argv[1]; tag = f[4:-4]
d = np.loadtxt(f)
t, Q, E, QS, R = d[:, 1], d[:, 4], d[:, 7], d[:, 10], d[:, 11]
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
ax[0].plot(t, Q, label="Q"); ax[0].plot(t, QS, "--", label="Q_S (sin zona fundida)")
for q in (-3, -2, -1, 0, 1, 2, 3): ax[0].axhline(q, color="0.88", lw=0.6)
ax[0].set_xlabel("t"); ax[0].set_ylabel("carga de Hopf"); ax[0].legend(fontsize=8)
ax[1].plot(t, R); ax[1].set_xlabel("t"); ax[1].set_ylabel("R_rms (peso 1 − n_z)")
ax[2].plot(t, E); ax[2].set_xlabel("t"); ax[2].set_ylabel("E2 + E4")
fig.suptitle(open(f).readline()[2:120], fontsize=7); fig.tight_layout()
fig.savefig(f"persistente_{tag}.png", dpi=110); print(f"persistente_{tag}.png")
