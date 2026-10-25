#!/usr/bin/env python3
# cola_angular.py — forma angular de la cola del hopfión  († Claude, 25-10-2026)
#
# Uso:  python3 cola_angular.py n_eta_N150_K1.bin  [otros .bin ...]
#
# Qué comprueba:
#   Lejos del núcleo el director está casi en ẑ y las desviaciones w = n_x + i·n_y cumplen
#   la ecuación lineal (∇² − k²) w = 0. Si sobrevive sólo el término dipolar,
#       w ∝ (x + i y)·k₁(kr)   →   |w|² ∝ sin²θ · (cola radial)
#   con θ medido desde el eje del nudo. Es decir: cero en el eje, máximo en el ecuador.
#
# Método:
#   1. Eje del nudo: autovector "distinto" del tensor de inercia pesado con (1 − n_z)
#      (un toro tiene dos momentos iguales y uno distinto). Así no se supone que el nudo
#      siga alineado con z.
#   2. Capas esféricas entre R_in (2,5 por defecto; RIN=3 python3 ... para cambiarlo) y 0,35·L. En cada capa se promedia
#      |w|² = n_x² + n_y² en bandas de cos θ y se normaliza por la media de la capa.
#   3. Se descompone en polinomios de Legendre:  f(μ) = a0 + a2·P2(μ) + a4·P4(μ) (+ impares).
#      sin²θ puro da  a2/a0 = −1  y  a4/a0 = 0.   Isótropo da a2/a0 = 0.
#   4. Control del régimen lineal: compara |w|² con 2(1 − n_z); deben coincidir si |w| ≪ 1.
#
# Salida: tabla por capa y la figura cola_angular.png.
# Necesita numpy y matplotlib.
import sys, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

def leer(fich):
    with open(fich, 'rb') as f:
        N = int(np.frombuffer(f.read(4), np.int32)[0]); L = float(np.frombuffer(f.read(8))[0])
        n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3)
    return N, L, n

def P(l, m):
    return {0: np.ones_like(m), 1: m, 2: 0.5*(3*m**2-1), 3: 0.5*(5*m**3-3*m),
            4: (35*m**4-30*m**2+3)/8}[l]

fichs = sys.argv[1:]
if not fichs:
    print(__doc__ if __doc__ else 'Uso: python3 cola_angular.py fichero.bin'); sys.exit(1)

fig, ax = plt.subplots(len(fichs), 2, figsize=(12, 4.6*len(fichs)), squeeze=False)
for fi, fich in enumerate(fichs):
    N, L, n = leer(fich)
    h = L/(N-1); x = -L/2 + h*np.arange(N)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r = np.sqrt(X**2+Y**2+Z**2)
    nz = n[..., 2]; w2 = n[..., 0]**2 + n[..., 1]**2; d = 1 - nz

    # 1. eje del nudo
    pos = np.stack([X, Y, Z], -1)
    P3 = pos.reshape(-1, 3); dd = d.ravel()
    I = (P3 * dd[:, None]).T @ P3 / dd.sum()
    ev, evec = np.linalg.eigh(I)
    dif = [abs(ev[i] - np.mean(np.delete(ev, i))) for i in range(3)]
    ieje = int(np.argmax(dif)); eje = evec[:, ieje]
    if eje[2] < 0: eje = -eje
    incl = np.degrees(np.arccos(min(1, abs(eje[2]))))
    # radio del núcleo (aprox.): radio rms pesado con d
    R_rms = np.sqrt((r**2*d).sum()/d.sum())
    mu = np.where(r > 0, (pos @ eje)/np.maximum(r, 1e-12), 0)

    nom = os.path.basename(fich)
    print(f"\n{nom}:  N={N}  L={L}  h={h:.4f}")
    print(f"  eje del nudo = ({eje[0]:+.3f}, {eje[1]:+.3f}, {eje[2]:+.3f})   inclinación respecto a z = {incl:.1f}°")
    print(f"  momentos de inercia: {ev[0]:.3f} {ev[1]:.3f} {ev[2]:.3f}   R_rms = {R_rms:.3f}")

    # 2-3. capas
    # R_rms queda inflado por la propia cola: no sirve para fijar la ventana.
    # Inicio fijo (núcleo del anillo ≈ 1,2 → 1,2 + 1,3), ajustable con RIN=...; final 0,35·L.
    R_in = float(os.environ.get('RIN', 2.5)); R_out = 0.35*L
    if R_in >= R_out - 2*h:
        print(f"  ventana vacía: R_in = {R_in} ≥ 0,35·L = {R_out:.2f}. Baja RIN o usa una caja mayor."); continue
    bordes = np.arange(R_in, R_out + 1e-9, max(2*h, (R_out-R_in)/8))
    nb = 24; mb = np.linspace(-1, 1, nb+1); mc = 0.5*(mb[1:]+mb[:-1])
    print(f"  {'r':>6s} {'<|w|²>':>11s} {'lineal':>7s} {'a1/a0':>7s} {'a2/a0':>7s} {'a3/a0':>7s} {'a4/a0':>7s}   (sin²θ puro: a2/a0 = −1, resto 0)")
    rs, a2s, a4s = [], [], []
    cmap = plt.get_cmap('viridis')
    for j in range(len(bordes)-1):
        sel = (r >= bordes[j]) & (r < bordes[j+1])
        if sel.sum() < 200: continue
        m_s, w_s, d_s = mu[sel], w2[sel], d[sel]
        media = w_s.mean()
        if media < 1e-28: continue
        lin = media / (2*d_s.mean())          # ≈ 1 en régimen lineal
        idx = np.clip(np.digitize(m_s, mb)-1, 0, nb-1)
        prof = np.array([w_s[idx == b].mean() if np.any(idx == b) else np.nan for b in range(nb)]) / media
        ok = ~np.isnan(prof)
        # Legendre por mínimos cuadrados sobre las bandas (bandas de igual Δμ = igual área en la esfera)
        A = np.stack([P(l, mc[ok]) for l in range(5)], 1)
        a = np.linalg.lstsq(A, prof[ok], rcond=None)[0]
        rm = 0.5*(bordes[j]+bordes[j+1])
        print(f"  {rm:6.2f} {media:11.3e} {lin:7.3f} {a[1]/a[0]:+7.3f} {a[2]/a[0]:+7.3f} {a[3]/a[0]:+7.3f} {a[4]/a[0]:+7.3f}")
        rs.append(rm); a2s.append(a[2]/a[0]); a4s.append(a[4]/a[0])
        ax[fi, 0].plot(mc[ok], prof[ok], 'o-', ms=3, color=cmap(j/max(1, len(bordes)-2)), label=f"r = {rm:.2f}")
    ax[fi, 0].plot(mc, 1.5*(1-mc**2), 'k--', lw=2, label='1,5·sin²θ (dipolo)')
    ax[fi, 0].set_xlabel('cos θ  (θ desde el eje del nudo)'); ax[fi, 0].set_ylabel('|w|² / media de la capa')
    ax[fi, 0].set_title(f'{nom}: forma angular por capas'); ax[fi, 0].legend(fontsize=7); ax[fi, 0].grid(alpha=.3)
    ax[fi, 1].plot(rs, a2s, 'o-', label='a2/a0'); ax[fi, 1].plot(rs, a4s, 's-', label='a4/a0')
    ax[fi, 1].axhline(-1, color='k', ls='--', lw=1, label='sin²θ: a2/a0 = −1')
    ax[fi, 1].axhline(0, color='gray', lw=.8)
    ax[fi, 1].set_xlabel('r'); ax[fi, 1].set_title('Coeficientes de Legendre frente a r')
    ax[fi, 1].legend(fontsize=8); ax[fi, 1].grid(alpha=.3)
plt.tight_layout(); plt.savefig('cola_angular.png', dpi=120); print('\n-> cola_angular.png')
