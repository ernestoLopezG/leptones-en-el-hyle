#!/usr/bin/env python3
# medir_tamano.py — tamaño del hopfión y forma de la cola, a partir de los .bin  († Claude, 21-10-2026)
# Uso:  python3 medir_tamano.py fich1.bin fich2.bin ...      (sin argumentos: todos los n_*.bin de la carpeta)
# Lee el formato de hopfion.cpp / hopfion_eta.cpp (int N, double L, N³ vectores n).
# Mide, en el plano ecuatorial z = 0 (promediando en el ángulo):
#   R_anillo : radio donde n_z es mínimo (núcleo del anillo), interpolado
#   R_rms    : radio cuadrático medio pesado con (1 − n_z)  (tamaño global, menos sensible a la malla)
#   cola     : fuera del anillo, ajusta ln(1 − n_z) frente a r (exponencial, da λ) y frente a ln r (potencia, da p)
#              y dice cuál ajusta mejor.  k_ef: ajuste de ln[r²(1−n_z)] = cte − 2k·r (cola de Yukawa); debe dar ≈ k.
# Ventana: desde R_anillo + 1,5 hasta 0,35·L (lejos del borde fijo).
# Necesita numpy y matplotlib (Debian: sudo apt install python3-numpy python3-matplotlib). Dibuja tamano_hopfion.png
import sys, glob, os, re, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

fichs = sys.argv[1:] or sorted(glob.glob('n_*.bin'))
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
print(f"{'fichero':28s} {'N':>4s} {'L':>5s} {'R_anillo':>9s} {'R_rms':>7s} {'lambda':>7s} {'p':>6s} {'k_ef':>6s}  mejor ajuste")
for fich in fichs:
    with open(fich, 'rb') as f:
        N = int(np.frombuffer(f.read(4), np.int32)[0]); L = float(np.frombuffer(f.read(8))[0])
        n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3)
    h = L/(N-1); x = -L/2 + h*np.arange(N)
    nz = n[..., 2]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij'); r3 = np.sqrt(X**2+Y**2+Z**2)
    w = 1 - nz
    R_rms = np.sqrt((r3**2*w).sum()/w.sum())
    # plano ecuatorial (índice central en el eje z; con N par se promedian los dos centrales)
    k0 = N//2
    plano = nz[:, :, k0] if N % 2 else 0.5*(nz[:, :, k0-1] + nz[:, :, k0])
    X2, Y2 = np.meshgrid(x, x, indexing='ij'); rho = np.sqrt(X2**2+Y2**2)
    bins = np.arange(0, L/2, h); idx = np.digitize(rho.ravel(), bins)
    prof = np.array([plano.ravel()[idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, len(bins))])
    rc = 0.5*(bins[1:]+bins[:-1]); ok = ~np.isnan(prof); rc, prof = rc[ok], prof[ok]
    i = np.argmin(prof)
    if 0 < i < len(prof)-1:   # parábola por tres puntos
        a, b, c = np.polyfit(rc[i-1:i+2], prof[i-1:i+2], 2); R_an = -b/(2*a)
    else: R_an = rc[i]
    sel = (rc > R_an + 1.5) & (rc < 0.35*L) & (1 - prof > 1e-12)
    r_t, y_t = rc[sel], np.log(1 - prof[sel])
    if len(r_t) > 4:
        pe = np.polyfit(r_t, y_t, 1); re_ = y_t - np.polyval(pe, r_t)
        pp = np.polyfit(np.log(r_t), y_t, 1); rp = y_t - np.polyval(pp, np.log(r_t))
        lam, p = -1/pe[0], -pp[0]
        # predicción de la cola linealizada con giro k:  1 − n_z ≈ C·e^(−2kr)/r²  →  ln[r²(1−n_z)] = cte − 2k·r
        py = np.polyfit(r_t, y_t + 2*np.log(r_t), 1); k_ef = -py[0]/2
        mejor = 'exponencial' if (re_**2).sum() < (rp**2).sum() else 'potencia'
        resid = f"(res. exp {np.sqrt((re_**2).mean()):.3f}, pot {np.sqrt((rp**2).mean()):.3f})"
    else: lam = p = k_ef = np.nan; mejor = 'ventana demasiado corta'; resid = ''
    nom = os.path.basename(fich)
    print(f"{nom:28s} {N:4d} {L:5.1f} {R_an:9.3f} {R_rms:7.3f} {lam:7.3f} {p:6.2f} {k_ef:6.3f}  {mejor} {resid}")
    m = re.search(r'K([\d.]+)', nom); et = f"k = {m.group(1)}" if m else 'k = 0 (hopfion.cpp)'
    ax[0].plot(rc, prof, label=f"{et}  ({nom})")
    ax[1].semilogy(rc, np.clip(1 - prof, 1e-14, None), label=et)
ax[0].set_xlabel('r en el plano z = 0'); ax[0].set_ylabel('n_z (promedio en el ángulo)')
ax[0].set_title('Perfil del director: el mínimo es el núcleo del anillo'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[1].set_xlabel('r'); ax[1].set_ylabel('1 − n_z  (escala log)')
ax[1].set_title('Cola: recta aquí = decaimiento exponencial'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig('tamano_hopfion.png', dpi=120); print('-> tamano_hopfion.png')
