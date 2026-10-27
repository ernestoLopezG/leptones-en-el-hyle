#!/usr/bin/env python3
# figura_espin_flujo.py — líneas de flujo del hopfión para espín +½, −½ y la imagen especular  († Claude, 26-10-2026)
#
# Campo analítico: mapa de Hopf de carga 1 (n → ẑ en el infinito, n = −ẑ en el anillo r = 1, z = 0).
#   W = z1/z2,  z1 = 2(x+iy)/(1+r²),  z2 = (2z + i(r²−1))/(1+r²),   n = (2 Re W, 2 Im W, 1−|W|²)/(1+|W|²)
# Vorticidad del flujo (relación de Mermin–Ho):  ω_i = ½ ε_ijk n·(∂_j n × ∂_k n), calculada por diferencias.
# Las líneas de ω son cerradas y se enlazan; el enlace se mide con la integral de Gauss.
#   +½ : el nudo tal cual.      −½ : el mismo nudo girado 180° alrededor de x  (n'(x) = n(R⁻¹x)).
#   espejo : reflejado en z → −z (n'(x) = n(Px)): no es un giro, cambia el signo del enlace.
# Uso: python3 figura_espin_flujo.py   → espin_flujo.png
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def n_hopf(p):
    x, y, z = p; r2 = x*x+y*y+z*z
    z1 = 2*(x+1j*y)/(1+r2); z2 = (2*z + 1j*(r2-1))/(1+r2)
    if abs(z2) < 1e-14: return np.array([0., 0., -1.])
    W = z1/z2; a = abs(W)**2
    return np.array([2*W.real, 2*W.imag, 1-a])/(1+a)

R = np.diag([1., -1., -1.])      # giro de 180° alrededor de x (R⁻¹ = R)
P = np.diag([1., 1., -1.])       # espejo z → −z
casos = {'+½': lambda p: n_hopf(p),
         '−½': lambda p: n_hopf(R @ p),
         'espejo': lambda p: n_hopf(P @ p)}

def vort(nf, p, h=1e-4):
    d = [(nf(p+h*e)-nf(p-h*e))/(2*h) for e in np.eye(3)]
    n0 = nf(p)
    return np.array([n0 @ np.cross(d[1], d[2]), n0 @ np.cross(d[2], d[0]), n0 @ np.cross(d[0], d[1])])

def linea(nf, p0, L=14.0):
    f = lambda s, p: (lambda w: w/np.linalg.norm(w))(vort(nf, p))
    def vuelta(s, p): return np.linalg.norm(p-p0) - 0.02
    vuelta.terminal = True
    sol = solve_ivp(f, [0, L], p0, max_step=0.03, rtol=1e-8, atol=1e-10, events=None, dense_output=False)
    pts = sol.y.T
    # cortar al cerrar la curva (primera vuelta)
    dist = np.linalg.norm(pts-p0, axis=1); s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
    k = np.where((s > 1.0) & (dist < 0.05))[0]
    return pts[:k[0]+1] if len(k) else pts

def gauss_link(a, b):
    da = np.diff(a, axis=0); db = np.diff(b, axis=0); am = 0.5*(a[1:]+a[:-1]); bm = 0.5*(b[1:]+b[:-1])
    tot = 0.0
    for i in range(len(am)):
        rr = am[i]-bm; d3 = np.linalg.norm(rr, axis=1)**3
        tot += np.sum(np.einsum('ij,ij->i', np.cross(da[i], db), rr)/d3)
    return tot/(4*np.pi)

semillas = [0.45, 0.7, 1.45, 2.2]      # puntos (ρ, 0, 0) en el plano del anillo; ρ = 1 es el núcleo
colores = ['#3d5a80', '#5b8e7d', '#e07a5f', '#b5838d']
fig = plt.figure(figsize=(17, 7.2))
titulos = {'+½': 'Espín +½', '−½': 'Espín −½\n(el mismo nudo girado 180° alrededor de x)',
           'espejo': 'Imagen especular de +½\n(no es un cambio de espín)'}
for j, (nom, nf) in enumerate(casos.items()):
    ax = fig.add_subplot(1, 3, j+1, projection='3d')
    curvas = []
    M = R if nom == '−½' else (P if nom == 'espejo' else np.eye(3))
    for s, c in zip(semillas, colores):
        p0 = M @ np.array([s, 0., 0.]) if nom == '−½' else np.array([s, 0., 0.])
        p0 = (R @ np.array([s, 0, 0.])) if nom == '−½' else ((P @ np.array([s, 0, 0.])) if nom == 'espejo' else np.array([s, 0, 0.]))
        cv = linea(nf, p0); curvas.append(cv)
        ax.plot(cv[:, 0], cv[:, 1], cv[:, 2], color=c, lw=2)
        for q in np.linspace(0, len(cv)-2, 4, dtype=int)[1:]:          # flechas de sentido
            d = cv[q+1]-cv[q]; d = d/np.linalg.norm(d)*0.35
            ax.quiver(cv[q, 0], cv[q, 1], cv[q, 2], d[0], d[1], d[2], color=c, arrow_length_ratio=0.6, lw=2)
    # núcleo del anillo (n = −ẑ): circunferencia r = 1 en el plano del nudo
    t = np.linspace(0, 2*np.pi, 200); ring = np.stack([np.cos(t), np.sin(t), 0*t], 1) @ M.T
    ax.plot(ring[:, 0], ring[:, 1], ring[:, 2], 'k--', lw=1)
    # chorro central: sentido de ω sobre el eje
    w0 = vort(nf, np.array([0, 0, 0.0])); w0 = w0/np.linalg.norm(w0)
    ax.quiver(0, 0, -1.8*w0[2], 0, 0, 3.6*w0[2], color='crimson', lw=3, arrow_length_ratio=0.12)
    ax.text(0.3, 0, 1.3*np.sign(w0[2]), 'eje y espín' if nom != 'espejo' else 'eje', color='crimson', fontsize=9)
    lk = gauss_link(curvas[0], curvas[2])
    ax.set_title(f"{titulos[nom]}\nenlace de dos líneas = {lk:+.2f}", fontsize=11, y=0.98)
    ax.set_xlim(-2.0, 2.0); ax.set_ylim(-2.0, 2.0); ax.set_zlim(-2.0, 2.0); ax.set_box_aspect((1, 1, 1), zoom=1.45)
    ax.view_init(elev=22, azim=-55); ax.set_axis_off()
    print(nom, 'enlace', round(lk, 3), 'ω en el centro', np.round(w0, 3))
plt.suptitle('Líneas de flujo (vorticidad de Mermin–Ho) del hopfión: girar el nudo invierte todas las líneas vistas desde fuera,\n'
             'pero conserva el enlace; solo el espejo lo cambia de signo (mapa de Hopf analítico, ilustrativo)', fontsize=11)
plt.subplots_adjust(left=0, right=1, bottom=0, top=0.80, wspace=0); plt.savefig('espin_flujo.png', dpi=130); print('-> espin_flujo.png')
