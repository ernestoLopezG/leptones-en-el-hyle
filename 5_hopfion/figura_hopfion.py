# figura_hopfion.py — figura del hopfión con la SOLUCIÓN NUMÉRICA (no esquema)  († Claude, 21-10-2026)
#   (a) superficies de igual n_z (toros anidados) en 3D, cortadas un cuarto, con líneas del flujo de Mermin–Ho
#       integradas desde puntos sobre esas superficies;
#   (b) corte por el plano y = 0: n_z en color, contornos de n_z y vectores del flujo en el plano.
# Uso: python3 figura_hopfion.py [n_final_N150.bin]      (sin argumento: el n_final_N*.bin más reciente de la carpeta)
# Requiere numpy, scipy, matplotlib y scikit-image (pip install scikit-image).
import numpy as np, sys, glob, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from skimage import measure
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp
if len(sys.argv) > 1: fich = sys.argv[1]
else:
    c = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'n_final_N*.bin')), key=os.path.getmtime)
    if not c: sys.exit("No encuentro n_final_N*.bin; pásalo como argumento.")
    fich = c[-1]
print("leyendo", fich)
f = open(fich, 'rb'); N = np.frombuffer(f.read(4), np.int32)[0]; L = np.frombuffer(f.read(8))[0]
n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3); h = L/(N-1); x = -L/2 + h*np.arange(N)
# flujo de Mermin–Ho: ω = B = ½ε n·(∂n×∂n), v = A (gauge de Coulomb, por FFT)
d = [np.gradient(n, h, axis=a) for a in range(3)]
F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
B = np.stack([F(1, 2), F(2, 0), F(0, 1)], -1)
k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij'); K2 = KX**2+KY**2+KZ**2; K2[0, 0, 0] = 1
Bk = [np.fft.fftn(B[..., i]) for i in range(3)]
v = np.stack([np.real(np.fft.ifftn(c)) for c in [1j*(KY*Bk[2]-KZ*Bk[1])/K2, 1j*(KZ*Bk[0]-KX*Bk[2])/K2, 1j*(KX*Bk[1]-KY*Bk[0])/K2]], -1)
nz = n[..., 2]
# centro y eje: el mínimo de n_z marca el núcleo del anillo
i0 = np.unravel_index(np.argmin(nz), nz.shape); print("mínimo de n_z =", nz[i0], "en", x[list(i0)])
fig = plt.figure(figsize=(15, 7))
ax = fig.add_subplot(1, 2, 1, projection='3d')
niveles = [(-0.6, '#c0392b'), (0.0, '#e67e22'), (0.6, '#2e86c1')]
s = slice(None)
for lev, col in niveles:
    verts, faces, _, _ = measure.marching_cubes(nz, lev, spacing=(h, h, h))
    verts = verts - L/2
    cen = verts[faces].mean(1); keep = ~((cen[:, 0] > 0) & (cen[:, 1] < 0))   # quita un cuarto
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces[keep], verts[:, 2], color=col, alpha=0.22, linewidth=0, shade=True)
# líneas de VORTICIDAD (ω = B): yacen sobre los toros y son las que se enlazan (helicidad = carga de Hopf)
vi = [RegularGridInterpolator((x, x, x), B[..., c], bounds_error=False, fill_value=0) for c in range(3)]
def campo(t, p): return [vi[c](p)[0] for c in range(3)]
iz = np.argmin(abs(x)); xs = x[iz:]; prof = nz[iz:, iz, iz]
for lev, col in niveles:
    j = np.argmin(abs(prof - lev)); r0 = xs[j]
    for ang in (0.0, 2.1):
        p0 = [r0*np.cos(ang), r0*np.sin(ang), 0.0]
        def vuelta(t, p): return np.hypot(p[0]-p0[0], np.hypot(p[1]-p0[1], p[2]-p0[2])) - 0.05 if t > 0.5 else 1.0
        vuelta.terminal = True
        sol = solve_ivp(campo, (0, 60), p0, max_step=0.02, rtol=1e-7, events=vuelta)
        P = sol.y; m = ~((P[0] > 0) & (P[1] < 0))
        P = np.where(m, P, np.nan)
        ax.plot(P[0], P[1], P[2], color=col, lw=1.3)
ax.plot([0, 0], [0, 0], [-2.5, 2.5], color='gray', ls='--', lw=1)
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_zlim(-2.5, 2.5); ax.set_box_aspect((1, 1, 0.85)); ax.view_init(28, -35); ax.set_axis_off()
ax.set_title('Solución numérica: superficies n_z = −0,6 / 0 / 0,6\n(toros anidados, cortados un cuarto) y líneas de vorticidad\ndel flujo de Mermin–Ho: circunferencias enlazadas', fontsize=10)
bx = fig.add_subplot(1, 2, 2)
j = N//2; ext = [x[0], x[-1], x[0], x[-1]]
im = bx.imshow(nz[:, j, :].T, origin='lower', extent=ext, cmap='RdBu', vmin=-1, vmax=1)
cs = bx.contour(x, x, nz[:, j, :].T, levels=[l for l, _ in niveles], colors=['k', 'k', 'k'], linewidths=[2, 1.5, 1], linestyles=['solid', 'dashed', 'dotted'])
bx.clabel(cs, fmt={-0.6: '−0,6', 0.0: '0', 0.6: '0,6'}, fontsize=9)
st = max(1, N//28); X, Z = np.meshgrid(x, x, indexing='ij')
bx.quiver(X[::st, ::st], Z[::st, ::st], v[::st, j, ::st, 0], v[::st, j, ::st, 2], color='k', width=0.003)
plt.colorbar(im, ax=bx, label='n_z (componente del director a lo largo del eje)')
bx.set_xlim(-3.5, 3.5); bx.set_ylim(-3.5, 3.5); bx.set_aspect('equal'); bx.set_xlabel('x'); bx.set_ylabel('z')
bx.set_title('Corte y = 0: director (color), toros como curvas cerradas\n(niveles −0,6 / 0 / 0,6) y flujo en el plano (flechas)', fontsize=10)
plt.tight_layout(); out = os.path.join(os.path.dirname(os.path.abspath(fich)), 'figura_hopfion.png'); plt.savefig(out, dpi=130)
print("figura guardada en", out)
