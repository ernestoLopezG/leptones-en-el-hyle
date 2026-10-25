#!/usr/bin/env python3
# figura_cola_dipolar.py — intuición visual de la cola de Yukawa dipolar  († Claude, 25-10-2026)
# Campo analítico (no la simulación): fuera del núcleo, n_⊥ = n_x + i n_y = A·(x+iy)/r · k₁(kr),
# con k₁(z) ∝ e^(−z)(1+1/z)/z²  →  |n_⊥|² ∝ sen²θ · e^(−2kr)(1+1/kr)²/r².
# Unidades: r en 1/k (la longitud de Compton si k = 1/ƛ_C). El nudo se dibuja muy pequeño (ancla).
# Uso: python3 figura_cola_dipolar.py   → cola_dipolar.png
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from skimage.measure import marching_cubes

k = 1.0; A = 0.35; Rnudo = 0.12           # radio del nudo (ancla), en unidades de 1/k
def g(r):                                   # parte radial: k₁(z) ∝ e^(−z)(1+1/z)/z, con z = kr
    z = k*r; return np.exp(-z)*(1+1/z)/z
def wmod2(x, y, z):
    r = np.sqrt(x*x+y*y+z*z) + 1e-9
    s2 = (x*x+y*y)/r**2
    return (A*g(r))**2 * s2

fig = plt.figure(figsize=(16, 7))
gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.2, 0.9])

# ---------- 1. superficies |n_⊥|² = cte en 3D, cortadas un cuarto, con el nudo en el centro
ax = fig.add_subplot(gs[:, 0], projection='3d')
M = 3.2; N = 121; xs = np.linspace(-M, M, N); X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')
F = wmod2(X, Y, Z)
corte = (X > 0) & (Y < 0)                                   # quitamos un cuarto para ver dentro
ls = LightSource(azdeg=315, altdeg=45)
colores = ['#f6bd60', '#e07a5f', '#3d5a80']; niveles = [2e-2, 2e-3, 2e-4]; alfas = [0.45, 0.35, 0.2]
for niv, col, al in zip(niveles, colores, alfas):
    G = np.where(corte, 0, F)
    v, f, *_ = marching_cubes(G, niv, step_size=1)
    v = v*(xs[1]-xs[0]) - M
    ax.plot_trisurf(v[:, 0], v[:, 1], f, v[:, 2], color=col, alpha=al, linewidth=0, shade=True, lightsource=ls, antialiased=True)
# nudo: pequeño toro
u, t = np.meshgrid(np.linspace(0, 2*np.pi, 40), np.linspace(0, 2*np.pi, 20))
a = 0.09; Rn3 = 0.3   # en 3D, tamaño exagerado para que se vea
ax.plot_surface((Rn3+a*np.cos(t))*np.cos(u), (Rn3+a*np.cos(t))*np.sin(u), a*np.sin(t), color='crimson', linewidth=0, zorder=10)
ax.plot([0, 0], [0, 0], [-M, M], 'k--', lw=1)
ax.text(0, 0, M*1.02, 'eje del nudo', fontsize=9, ha='center')
ax.set_xlim(-M, M); ax.set_ylim(-M, M); ax.set_zlim(-M, M); ax.set_box_aspect((1, 1, 1))
ax.view_init(elev=22, azim=-60); ax.dist = 6.5; ax.set_axis_off()
ax.set_title('1. Superficies de igual deformación |n⊥|²\n(cortadas un cuarto; en rojo, el nudo, ampliado)', fontsize=10)

# ---------- 2. corte meridiano: color = |n_⊥|² (log), flechas = director en el plano xz
ax2 = fig.add_subplot(gs[:, 1])
M2 = 3.2; x2 = np.linspace(-M2, M2, 401); XX, ZZ = np.meshgrid(x2, x2, indexing='xy')
F2 = wmod2(XX, 0*XX, ZZ)
im = ax2.pcolormesh(XX, ZZ, np.log10(np.clip(F2, 1e-6, None)), cmap='magma', shading='auto', vmin=-6, vmax=-0.5)
plt.colorbar(im, ax=ax2, fraction=0.046, label='log₁₀ |n⊥|²')
# flechas del director (en y = 0, n_⊥ es real: n_x)
xa = np.linspace(-3, 3, 13); XA, ZA = np.meshgrid(xa, xa, indexing='xy')
RA = np.sqrt(XA**2+ZA**2)+1e-9
nx = 6*A*g(RA)*XA/RA; nx = np.clip(nx, -0.8, 0.8); nz = np.sqrt(1-nx**2)   # inclinación ×6 para que se vea
m = RA > 0.45
ax2.quiver(XA[m], ZA[m], nx[m], nz[m], color='white', scale=18, width=0.004, headwidth=3, pivot='middle')
for c in (1, 2, 3):
    ax2.add_patch(plt.Circle((0, 0), c, fill=False, ls=':', color='w', lw=0.8))
    ax2.text(c*0.71, c*0.71, f'{c}/k', color='w', fontsize=8)
ax2.plot([-Rnudo, Rnudo], [0, 0], 'o', color='crimson', ms=4)
ax2.axvline(0, color='w', ls='--', lw=0.8)
ax2.set_aspect('equal'); ax2.set_xlabel('x  (en unidades de 1/k)'); ax2.set_ylabel('z (eje del nudo)')
ax2.set_title('2. Corte por el eje\ncolor: deformación; flechas: director (inclinación ×6)', fontsize=10)

# ---------- 3. diagrama polar: sen²θ, y comparación de la caída radial
ax3 = fig.add_subplot(gs[0, 2], projection='polar')
th = np.linspace(0, 2*np.pi, 400)
ax3.plot(th, np.sin(th)**2, color='#3d5a80', lw=2)
ax3.set_theta_zero_location('N'); ax3.set_rticks([]); ax3.set_xticks(np.radians([0, 90, 180, 270]))
ax3.set_xticklabels(['eje', 'ecuador', 'eje', 'ecuador'], fontsize=8)
ax3.set_title('3a. Forma angular: sen²θ\n(cero en el eje, máxima en el ecuador)', fontsize=10)
ax4 = fig.add_subplot(gs[1, 2])
r = np.linspace(0.3, 6, 300)
ax4.semilogy(r, (g(r))**2, color='#e07a5f', lw=2, label='halo del director: e^(−2kr)/r² (Yukawa)')
ax4.semilogy(r, 0.5/r**4, color='#3d5a80', lw=2, ls='--', label='campo externo: 1/r⁴ ∝ |E|² (Coulomb)')
ax4.axvline(1, color='gray', lw=0.8); ax4.text(1.05, 1e-6, '1/k', fontsize=8)
ax4.set_xlabel('r (en 1/k)'); ax4.set_ylim(1e-8, 1e2); ax4.legend(fontsize=7); ax4.grid(alpha=.3)
ax4.set_title('3b. Caída radial: el halo se apaga, el campo externo no', fontsize=10)

plt.suptitle('La cola de Yukawa dipolar alrededor de un nudo pequeño (campo analítico, ilustrativo)', fontsize=12)
plt.tight_layout(); plt.savefig('cola_dipolar.png', dpi=130); print('-> cola_dipolar.png')
