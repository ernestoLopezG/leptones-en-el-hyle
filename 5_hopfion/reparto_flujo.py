#!/usr/bin/env python3
# reparto_flujo.py — reparto de la energía del flujo de Mermin–Ho entre poloidal y toroidal  († Claude, 27-10-2026)
#
# Pregunta: ¿el reparto 66 % / 34 % del hopfión numérico es 2/3 y 1/3?
# Mide la fracción toroidal Et/E (la que gira alrededor del eje) en:
#   (a) el mapa de Hopf analítico (sin error de minimización), a varias mallas y cajas, para ver a qué converge;
#   (b) los .bin que se pasen como argumento (hopfión numérico de hopfion.cpp / hopfion_eta.cpp).
# El flujo se obtiene como en flujo.py: ω = B (Mermin–Ho), v = A en gauge de Coulomb (por FFT).
# Uso:  python3 reparto_flujo.py                      → solo el analítico
#       python3 reparto_flujo.py n_final_N100.bin n_final_N150.bin   → además los numéricos
import sys, numpy as np

def reparto(n, h):
    N = n.shape[0]; x = -(N-1)*h/2 + h*np.arange(N)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    d = [np.gradient(n, h, axis=a) for a in range(3)]
    F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
    B = np.stack([F(1, 2), F(2, 0), F(0, 1)], -1)
    k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    K2 = KX**2+KY**2+KZ**2; K2[0, 0, 0] = 1
    Bk = [np.fft.fftn(B[..., i]) for i in range(3)]
    v = np.stack([np.real(np.fft.ifftn(c)) for c in
                  [1j*(KY*Bk[2]-KZ*Bk[1])/K2, 1j*(KZ*Bk[0]-KX*Bk[2])/K2, 1j*(KX*Bk[1]-KY*Bk[0])/K2]], -1)
    R = np.sqrt(X**2+Y**2) + 1e-12
    vphi = (-Y*v[..., 0] + X*v[..., 1])/R
    E = (v**2).sum(); Et = (vphi**2).sum()
    Q = (v*B).sum()*h**3/(16*np.pi**2)
    return Et/E, Q

def hopf_analitico(N, L, a=1.0):
    x = -L/2 + L/(N-1)*np.arange(N); X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    X, Y, Z = X/a, Y/a, Z/a; r2 = X**2+Y**2+Z**2
    z1 = 2*(X+1j*Y)/(1+r2); z2 = (2*Z + 1j*(r2-1))/(1+r2)
    W = z1/np.where(np.abs(z2) < 1e-14, 1e-14, z2); A2 = np.abs(W)**2
    n = np.stack([2*W.real, 2*W.imag, 1-A2], -1)/(1+A2)[..., None]
    return n, L/(N-1)

print("(a) mapa de Hopf analítico (radio del anillo a = 1)")
print(f"{'N':>5s} {'L':>5s} {'h':>7s} {'toroidal':>9s} {'poloidal':>9s} {'Q':>7s}")
for N, L in [(64, 12), (96, 12), (128, 12), (128, 20), (160, 20), (192, 24)]:
    n, h = hopf_analitico(N, L)
    ft, Q = reparto(n, h)
    print(f"{N:5d} {L:5.0f} {h:7.3f} {ft:9.4f} {1-ft:9.4f} {Q:7.3f}")

for fich in sys.argv[1:]:
    with open(fich, 'rb') as f:
        N = int(np.frombuffer(f.read(4), np.int32)[0]); L = float(np.frombuffer(f.read(8))[0])
        n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3)
    ft, Q = reparto(n, L/(N-1))
    print(f"(b) {fich}: N={N} L={L}  toroidal {ft:.4f}  poloidal {1-ft:.4f}  Q {Q:.3f}")
print("2/3 = 0.6667, 1/3 = 0.3333")
