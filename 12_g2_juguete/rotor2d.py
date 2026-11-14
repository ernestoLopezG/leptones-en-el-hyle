#!/usr/bin/env python3
"""
Juguete 2D para g-2: una fuente que gira en un medio donde las
perturbaciones viajan a c=1.  phi_tt = lap(phi) + s(x,t)
Se mide el momento angular del campo dentro de un radio r y se
compara con el que tendria el campo estatico arrastrado rigidamente
(cuasiestatico).  Cociente = fraccion del campo que acompaña a la fuente.
Uso: python3 rotor2d.py [v1 v2 ...]   (v = a*Omega, por defecto 0.1 0.2 0.4)
"""
import sys, numpy as np
a, sig, Q = 8.0, 1.5, 1.0          # radio de giro, anchura de la fuente, carga
N, dx, dt = 440, 1.0, 0.5
esp = 60                            # esponja absorbente
x = (np.arange(N)-N/2+0.5)*dx
X, Y = np.meshgrid(x, x, indexing='ij')
R = np.hypot(X, Y)
d = np.clip((R-(N/2-esp))/esp, 0, None); gam = 0.3*d**2

def fuente(t, Om, ramp):
    xs, ys = a*np.cos(Om*t), a*np.sin(Om*t)
    f = min(1.0, t/ramp)
    return f*Q/(2*np.pi*sig**2)*np.exp(-((X-xs)**2+(Y-ys)**2)/(2*sig**2)), xs, ys

def grad(p):
    gx = (np.roll(p,-1,0)-np.roll(p,1,0))/(2*dx)
    gy = (np.roll(p,-1,1)-np.roll(p,1,1))/(2*dx)
    return gx, gy

def lap(p):
    return (np.roll(p,1,0)+np.roll(p,-1,0)+np.roll(p,1,1)+np.roll(p,-1,1)-4*p)/dx**2

radios = np.array([10,15,20,30,40,60,80,100,130,150])

def corre(v):
    Om = v/a; T = 2*np.pi/Om
    ramp = 1.5*T; tfin = ramp + 300 + 2*T
    p0 = np.zeros((N,N)); p1 = np.zeros((N,N)); t = 0.0
    Lacum = np.zeros(len(radios)); Lqs = np.zeros(len(radios)); n = 0
    Pflux = 0.0
    while t < tfin:
        s, xs, ys = fuente(t, Om, ramp)
        p2 = (2*p1 - (1-gam*dt/2)*p0 + dt**2*(lap(p1)+s))/(1+gam*dt/2)
        t += dt
        if t > tfin - T:            # promedio sobre la ultima vuelta
            pt = (p2-p0)/(2*dt); gx, gy = grad(p1)
            ell = -pt*(X*gy - Y*gx)
            # campo estatico analitico de la fuente gaussiana, arrastrado
            xs, ys = a*np.cos(Om*(t-dt)), a*np.sin(Om*(t-dt))
            rx, ry = X-xs, Y-ys; rho = np.hypot(rx, ry)+1e-12
            Er = -Q*(1-np.exp(-rho**2/(2*sig**2)))/(2*np.pi*rho)   # dphi/drho
            sx, sy = Er*rx/rho, Er*ry/rho
            Vx, Vy = -v*np.sin(Om*(t-dt)), v*np.cos(Om*(t-dt))
            pts = -(Vx*sx+Vy*sy)
            ellqs = -pts*(X*sy - Y*sx)
            for i, r in enumerate(radios):
                m = R < r
                Lacum[i] += ell[m].sum()*dx**2; Lqs[i] += ellqs[m].sum()*dx**2
            # flujo de energia radiado a traves de r=150
            anillo = (R > 148) & (R < 152)
            Pflux += (-pt*(X*gx+Y*gy)/R)[anillo].sum()*dx**2/4.0  # /espesor
            n += 1
        p0, p1 = p1, p2
    return Om, Lacum/n, Lqs/n, Pflux/n

vs = [float(z) for z in sys.argv[1:]] or [0.1, 0.2, 0.4]
for v in vs:
    Om, L, Lq, P = corre(v)
    print(f"\nv = {v}   Omega = {Om:.4f}   c/Omega = {1/Om:.0f}   (a = {a}, sigma = {sig})")
    print("   r    r*Om   L_campo     L_cuasiest   fraccion")
    for r, l1, l2 in zip(radios, L, Lq):
        print(f"{r:5.0f} {r*Om:6.2f} {l1:11.4e} {l2:11.4e} {l1/l2:8.3f}")
