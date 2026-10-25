# flujo.py — Paso 2a: el flujo que genera el hopfión según la relación de Mermin–Ho  († Claude, 19-10-2026)
#   Hipótesis [H]: la vorticidad del flujo del hyle es proporcional al retorcimiento del director,
#       ω_i = ∇×v = ½ ε_ijk n·(∂_j n × ∂_k n) = B_i      (constante de proporcionalidad = 1 en estas unidades)
#   Entonces v = A (potencial vector en gauge de Coulomb, ∇·v = 0: flujo incompresible) y la
#   helicidad del flujo ∫ v·ω d³x = 16π² Q  (Q = carga de Hopf): el enlace del flujo ES la topología del director.
# Uso: python3 flujo.py n_final_N<N>.bin      -> flujo.png y cifras en pantalla
import numpy as np, sys, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import glob, os
if len(sys.argv) > 1: fich = sys.argv[1]
else:   # sin argumento: usa el n_final_N*.bin más reciente de la carpeta del programa
    cands = sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'n_final_N*.bin')), key=os.path.getmtime)
    if not cands: sys.exit("No encuentro ningún n_final_N*.bin en la carpeta de flujo.py; pásalo como argumento.")
    fich = cands[-1]
print("leyendo", fich)
f = open(fich, 'rb'); N = np.frombuffer(f.read(4), np.int32)[0]; L = np.frombuffer(f.read(8))[0]
n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3); h = L/(N-1)
x = -L/2 + h*np.arange(N); X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
d = [np.gradient(n, h, axis=a) for a in range(3)]
F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
B = np.stack([F(1, 2), F(2, 0), F(0, 1)], -1)                     # vorticidad ω
k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij'); K2 = KX**2+KY**2+KZ**2; K2[0,0,0] = 1
Bk = [np.fft.fftn(B[..., i]) for i in range(3)]
v = np.stack([np.real(np.fft.ifftn(c)) for c in [1j*(KY*Bk[2]-KZ*Bk[1])/K2, 1j*(KZ*Bk[0]-KX*Bk[2])/K2, 1j*(KX*Bk[1]-KY*Bk[0])/K2]], -1)
dV = h**3
Hel = (v*B).sum()*dV; Q = Hel/(16*np.pi**2)
R = np.sqrt(X**2+Y**2) + 1e-12
vphi = (-Y*v[...,0] + X*v[...,1])/R                               # componente que gira alrededor del eje z (toroidal)
wphi = (-Y*B[...,0] + X*B[...,1])/R
Ez = ((X*v[...,1] - Y*v[...,0])).sum()*dV                          # "momento angular" del flujo (densidad 1) alrededor de z
Ekin = 0.5*(v**2).sum()*dV
print(f"Q = {Q:.4f}   helicidad ∫v·ω = {Hel:.2f}   energía cinética del flujo ½∫v² = {Ekin:.3f}")
print(f"momento angular del flujo alrededor de z, ∫(x v_y − y v_x) = {Ez:.3f}")
# energía de la parte poloidal frente a la toroidal
Ep = 0.5*((v**2).sum(-1) - vphi**2).sum()*dV; Et = 0.5*(vphi**2).sum()*dV
print(f"reparto de la energía del flujo: toroidal (gira alrededor de z) {Et/Ekin:.2f}, poloidal (alrededor del anillo) {Ep/Ekin:.2f}")
j = N//2; s = slice(None, None, max(1, N//30))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
ext = [x[0], x[-1], x[0], x[-1]]
im0 = ax[0].imshow(n[:, j, :, 2].T, origin='lower', extent=ext, cmap='RdBu', vmin=-1, vmax=1)
ax[0].set_title('Director: componente n_z (corte y = 0)\nn_z = −1 en el núcleo del anillo', fontsize=9); plt.colorbar(im0, ax=ax[0])
im1 = ax[1].imshow(vphi[:, j, :].T, origin='lower', extent=ext, cmap='PuOr')
ax[1].quiver(X[s, j, s], Z[s, j, s], v[s, j, s, 0], v[s, j, s, 2], color='k', scale=None, width=0.004)
ax[1].set_title('Flujo de Mermin–Ho (corte y = 0)\nflechas: poloidal, alrededor del núcleo; color: toroidal, alrededor de z', fontsize=9); plt.colorbar(im1, ax=ax[1])
im2 = ax[2].imshow(np.sqrt((B[:, j, :]**2).sum(-1)).T, origin='lower', extent=ext, cmap='magma')
ax[2].set_title('|vorticidad| = retorcimiento del director', fontsize=9); plt.colorbar(im2, ax=ax[2])
for a in ax: a.set_xlim(-4, 4); a.set_ylim(-4, 4); a.set_xlabel('x'); a.set_ylabel('z')
plt.tight_layout(); salida = os.path.join(os.path.dirname(os.path.abspath(fich)), 'flujo.png'); plt.savefig(salida, dpi=110); print('gráfica guardada en', salida)
