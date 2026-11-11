"""Figura de ensamblaje: qué hay en el disco compacto en cada punto del espacio extendido.
Esquema, no a escala. Uso: python3 fig00_ensamblaje.py  ->  ../figuras/fig00_ensamblaje.png"""
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, ConnectionPatch

fig = plt.figure(figsize=(12, 6.2))
ax = fig.add_axes([0.0, 0.0, 0.56, 1.0]); ax.set_xlim(-5.2, 5.2); ax.set_ylim(-5.2, 5.2)
ax.set_aspect('equal'); ax.axis('off')
# envolvente (orbital s)
X, Y = np.meshgrid(np.linspace(-5.2, 5.2, 400), np.linspace(-5.2, 5.2, 400))
R = np.hypot(X, Y); dens = np.exp(-2*R/1.6)
ax.imshow(dens, extent=[-5.2, 5.2, -5.2, 5.2], origin='lower', cmap='Blues', alpha=0.85, vmin=0, vmax=1)
ax.add_patch(Circle((0, 0), 2.6, fill=False, ls='--', lw=0.8, color='steelblue'))
ax.text(1.95, 2.25, 'envolvente $\\Psi(\\mathbf{x})$\n(orbital)', color='navy', fontsize=10, bbox=dict(fc='white', ec='none', alpha=0.75, pad=1))
# campo más allá: circulación u ~ 1/r (flechas en eta, esquemáticas) y silbido (anillos ondulados)
verde = 'seagreen'
for r in [1.3, 2.1, 3.0, 3.45, 3.9, 4.35]:
    t = np.linspace(0, 2*np.pi, 600)
    rr = r + 0.045*np.sin(36*t)
    ax.plot(rr*np.cos(t), rr*np.sin(t), color=verde, lw=1.4*min(1, 3.1/r), alpha=(0.35 if r < 2.6 else min(1, 0.9*3.1/r)))
ax.text(-5.05, -5.0, 'campo: circulación secundaria con signo, potencial $u\\propto 1/r$, que silba\n(tono universal, en el corte: no radia; no se sostiene solo)', color=verde, fontsize=9.5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
ax.text(-5.05, 4.75, 'pulsación (la partícula): azul, rojo y naranja; nudo: morado', color='navy', fontsize=9.5, bbox=dict(fc='white', ec='none', alpha=0.85, pad=1.5))
# nudo (toro esquemático) desplazado del centro: el nudo recorre la envolvente
kx, ky = 0.55, -0.35
ax.add_patch(Ellipse((kx, ky), 0.9, 0.42, fill=False, lw=2.2, color='purple'))
ax.add_patch(Ellipse((kx, ky), 0.32, 0.13, fill=False, lw=1.2, color='purple'))
ax.plot([kx, kx], [ky-0.55, ky+0.55], color='k', lw=1.2)
ax.text(kx+0.55, ky-0.75, 'nudo: centro del colapso y eje\ndel espín ($\\sim\\sqrt{3}\\,\\xi_0$)', color='purple', fontsize=9, bbox=dict(fc='white', ec='none', alpha=0.75, pad=1))

# discos compactos en tres puntos
puntos = [((0.4, 1.1), 1.0, 'dentro del orbital'), ((1.7, -1.9), 0.45, 'borde del orbital'), ((4.3, 0.2), 0.0, 'fuera: solo campo')]
pos_ins = [[0.60, 0.62, 0.20, 0.34], [0.60, 0.14, 0.20, 0.34], [0.80, 0.38, 0.20, 0.34]]
for (p, amp, tit), rect in zip(puntos, pos_ins):
    ax.plot(*p, 'ko', ms=4)
    a2 = fig.add_axes(rect); a2.set_xlim(-1.3, 1.3); a2.set_ylim(-1.45, 1.3); a2.set_aspect('equal'); a2.axis('off')
    a2.add_patch(Circle((0, 0), 1.15, fill=False, lw=0.8))                 # R_u
    a2.add_patch(Circle((0, 0), 0.16, color='k'))                           # núcleo fluido
    a2.add_patch(Circle((0, 0), 0.42, fill=False, lw=0.8, color='gray'))    # interfaz A/C*
    if amp > 0:
        t = np.linspace(0, 2*np.pi, 300)
        rr = 0.42 + 0.14*amp*np.cos(t - 0.6)**2
        a2.plot(rr*np.cos(t), rr*np.sin(t), color='navy', lw=1.0+1.5*amp)
        for z in [0.6, 0.6+np.pi]:
            a2.plot(0.42*np.cos(z), 0.42*np.sin(z), 'o', color='darkorange', ms=4+3*amp)
        a2.add_patch(FancyArrowPatch((0.62*np.cos(1.2), 0.62*np.sin(1.2)), (0.62*np.cos(2.1), 0.62*np.sin(2.1)),
                     connectionstyle='arc3,rad=0.35', arrowstyle='->', mutation_scale=10, color='crimson'))
    else:
        a2.add_patch(FancyArrowPatch((0.95*np.cos(1.0), 0.95*np.sin(1.0)), (0.95*np.cos(2.2), 0.95*np.sin(2.2)),
                     connectionstyle='arc3,rad=0.3', arrowstyle='->', mutation_scale=10, color='seagreen'))
        t = np.linspace(0, 2*np.pi, 500); rr = 0.42 + 0.04*np.sin(24*t)
        a2.plot(rr*np.cos(t), rr*np.sin(t), color='seagreen', lw=1.2)
    a2.set_title(tit, fontsize=9, pad=1)
    con = ConnectionPatch(xyA=p, coordsA=ax.transData, xyB=(0, -1.15), coordsB=a2.transData, color='0.5', lw=0.6, ls=':')
    fig.add_artist(con)
fig.text(0.585, 0.025,
 'Disco compacto $(\\xi,\\eta)$ en cada punto: núcleo fluido (negro), interfaz A/C* (gris),\n'
 'circulación de orden ½ (flecha roja); la curva azul es la intensidad\n'
 'de la onda a lo largo de η, '
 'con máximos en las zonas de mayor índice (naranja),\n'
 'escalada con $|\\Psi(\\mathbf{x})|$. El campo (verde) existe también dentro\n'
 'del orbital; fuera de él solo hay campo: la circulación secundaria\n'
 'con signo (flecha) y su silbido (ondulado), en la interfaz.', fontsize=8.5)
plt.savefig('../figuras/fig00_ensamblaje.png', dpi=160)
print('-> ../figuras/fig00_ensamblaje.png')
