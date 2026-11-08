# figuras_nucleacion.py — gráficas de nucleacion.cpp  (Claude, 2-11-2026)
#
# Uso:
#   python3 figuras_nucleacion.py <tag>            → nucleacion_<tag>.png: series temporales y cortes y = 0
#   python3 figuras_nucleacion.py resumen [tags…]  → nucleacion_resumen.png: Q(t) de varias ejecuciones juntas
#                                                     (sin tags: todas las log_*.dat de la carpeta)
import sys, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def lee_log(tag):
    d = np.loadtxt(f"log_{tag}.dat", comments="#")
    cab = open(f"log_{tag}.dat").readline().strip("# \n")
    return d, cab

def fases(ax, d):
    t, f = d[:, 1], d[:, 9]
    for v in (1, 2):
        idx = np.where(f > v)[0]
        if len(idx): ax.axvline(t[idx[0]], color="0.7", lw=0.8, ls=":")

def figura(tag):
    d, cab = lee_log(tag)
    t, Tmax, lam, Q, E24, EB = d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 7], d[:, 8]
    cortes = sorted(glob.glob(f"corte_{tag}_*.dat"))
    elige = [cortes[i] for i in np.linspace(0, len(cortes)-1, min(4, len(cortes))).astype(int)] if cortes else []
    fig = plt.figure(figsize=(14, 8.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1.25])
    a = fig.add_subplot(gs[0, :2]); a.plot(t, Tmax, label="T máx (fusión si > 1)")
    a.plot(t, lam/max(lam.max(), 1e-9), label="λ (relativo)"); a.axhline(1, color="r", lw=0.6, ls="--")
    a.set_ylabel("protocolo"); a.legend(fontsize=8); fases(a, d)
    b = fig.add_subplot(gs[0, 2:]); b.plot(t, Q, color="k"); b.set_ylabel("Q (carga de Hopf)")
    for y in (-1, 0, 1): b.axhline(y, color="0.8", lw=0.6)
    b.set_title("Q solo es entero cuando todo el medio está ordenado (T = 0)", fontsize=8); fases(b, d)
    c = fig.add_subplot(gs[1, :2]); c.plot(t, E24, color="C2"); c.axhline(544, color="C2", ls="--", lw=0.7)
    c.set_ylabel("E2 + E4"); c.text(t[-1], 544, " hopfión ≈ 544", va="bottom", ha="right", fontsize=8); fases(c, d)
    c.set_xlabel("t")
    e = fig.add_subplot(gs[1, 2:]); e.plot(t, EB, color="C3"); e.set_ylabel("∫ S³ Ω·B (alineación con el flujo)")
    e.set_xlabel("t"); fases(e, d)
    for i, fn in enumerate(elige):
        tt = float(re.search(r"t=([0-9.]+)", open(fn).readline()).group(1))
        x, z, S, nz = np.loadtxt(fn, unpack=True); M = int(round(np.sqrt(len(x))))
        ax = fig.add_subplot(gs[2, i])
        im = ax.imshow(nz.reshape(M, M).T, origin="lower", extent=[x.min(), x.max(), z.min(), z.max()],
                       cmap="RdBu", vmin=-1, vmax=1)
        ax.contour(x.reshape(M, M), z.reshape(M, M), S.reshape(M, M), levels=[0.05], colors="k", linewidths=0.8)
        ax.set_title(f"t = {tt:.3f}: n_z en y = 0", fontsize=8); ax.set_xlabel("x"); ax.set_ylabel("z")
    if elige: fig.colorbar(im, ax=fig.axes[-1], fraction=0.046, label="n_z  (línea negra: borde de la zona fundida)")
    fig.suptitle(cab, fontsize=7); fig.tight_layout()
    fig.savefig(f"nucleacion_{tag}.png", dpi=110); print(f"guardada nucleacion_{tag}.png")

def resumen(tags):
    if not tags: tags = sorted(f[4:-4] for f in glob.glob("log_*.dat"))
    fig, ax = plt.subplots(figsize=(10, 5))
    for tg in tags:
        d, cab = lee_log(tg)
        ax.plot(d[:, 1], d[:, 4], label=tg, lw=1)
    for y in (-2, -1, 0, 1, 2): ax.axhline(y, color="0.85", lw=0.6)
    ax.set_xlabel("t"); ax.set_ylabel("Q"); ax.legend(fontsize=7, ncol=2)
    ax.set_title("Carga de Hopf durante la fusión, el enfriamiento y la retirada del flujo")
    fig.tight_layout(); fig.savefig("nucleacion_resumen.png", dpi=110); print("guardada nucleacion_resumen.png")

if __name__ == "__main__":
    if len(sys.argv) < 2: print(__doc__ if __doc__ else "uso: figuras_nucleacion.py <tag> | resumen [tags]"); sys.exit(1)
    if sys.argv[1] == "resumen": resumen(sys.argv[2:])
    else: figura(sys.argv[1])
