#!/usr/bin/env python3
# precesion_hopfion.py — hopfión entero (superficie n_z = 0 de la solución numérica) con su eje inclinado
# y precesionando alrededor de un campo externo  († Claude, 25-10-2026)
# Uso: python3 precesion_hopfion.py n_final_N150.bin [angulo_cono_grados=30]
# Necesita numpy, matplotlib y scikit-image (Debian: sudo apt install python3-skimage)
import sys, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes
fich = sys.argv[1]; th = np.radians(float(sys.argv[2]) if len(sys.argv) > 2 else 30)
with open(fich, 'rb') as f:
    N = int(np.frombuffer(f.read(4), np.int32)[0]); L = float(np.frombuffer(f.read(8))[0])
    n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3)
h = L/(N-1)
# superficie n_z = 0 (toro exterior), malla aligerada para que el sombreado salga limpio
v, fc, nor, _ = marching_cubes(n[..., 2], level=0.0, spacing=(h, h, h), step_size=2)
v -= L/2
# girar el toro: su eje (z) pasa a apuntar a lo largo de una generatriz del cono (ángulo th con z, en el plano xz)
R = np.array([[np.cos(th), 0, np.sin(th)], [0, 1, 0], [-np.sin(th), 0, np.cos(th)]])
v = v @ R.T; nor = nor @ R.T
# sombreado lambertiano a mano (evita el «ruido» del sombreado automático)
luz = np.array([0.4, -0.6, 0.7]); luz /= np.linalg.norm(luz)
nf = nor[fc].mean(axis=1); nf /= np.linalg.norm(nf, axis=1)[:, None]
I = 0.35 + 0.65*np.clip(np.abs(nf @ luz), 0, 1)
col = np.c_[0.12*I + 0.05, 0.45*I + 0.05, 0.75*I + 0.1, np.ones_like(I)]
fig = plt.figure(figsize=(8, 9)); ax = fig.add_subplot(projection='3d')
ax.add_collection3d(Poly3DCollection(v[fc], facecolors=col, edgecolor='none'))
# cono de precesión, eje actual y campo externo
H = 1.6*np.abs(v).max(); ph = np.linspace(0, 2*np.pi, 200); r = H*np.tan(th)
ax.plot(r*np.cos(ph), r*np.sin(ph), H*np.ones_like(ph), 'k', lw=2)
for p in np.linspace(0, 2*np.pi, 13)[:-1]:
    ax.plot([0, r*np.cos(p)], [0, r*np.sin(p)], [0, H], color='0.6', lw=0.6)
ax.plot([0, r], [0, 0], [0, H], 'k', lw=3)                                   # eje del espín (generatriz actual)
ax.quiver(0, 0, 0, 0, 0, 1.25*H, color='crimson', lw=2.5, arrow_length_ratio=0.08)   # campo externo
ax.text(0, 0, 1.3*H, 'campo externo', color='crimson', ha='center')
ax.text(1.15*r, 0, 0.8*H, 'eje del espín', ha='left')
a = np.linspace(0.15, 1.1, 40); ax.plot(r*np.cos(a), r*np.sin(a), H*np.ones_like(a)+0.02*H, color='crimson', lw=2)
ax.set_box_aspect((1, 1, 1)); m = 1.1*H
ax.set_xlim(-m/1.5, m/1.5); ax.set_ylim(-m/1.5, m/1.5); ax.set_zlim(-0.5*m, 1.2*m)
ax.view_init(elev=14, azim=-60); ax.set_axis_off()
ax.set_title('Precesión del eje del espín bajo un campo externo\n(hopfión numérico, superficie n_z = 0, entero)', fontsize=11)
plt.tight_layout(); plt.savefig('precesion_hopfion.png', dpi=130, bbox_inches='tight'); print('-> precesion_hopfion.png')
