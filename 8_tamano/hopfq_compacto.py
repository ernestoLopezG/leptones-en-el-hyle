#!/usr/bin/env python3
# hopfq_compacto.py — carga de Hopf Q = (1/16π²) ∫ A·B d³x de los .bin de hopfion_compacto.cpp  (Claude, 4-11-2026)
# z periódica (derivadas con np.roll); x, y con borde fijo (np.gradient). A en gauge de Coulomb por FFT.
# Uso:  python3 hopfq_compacto.py n_c_Lc*.bin
import sys, glob, numpy as np
for fich in sys.argv[1:] or sorted(glob.glob('n_*_Lc*.bin')):
    with open(fich, 'rb') as f:
        N = int(np.frombuffer(f.read(4), np.int32)[0]); Nz = int(np.frombuffer(f.read(4), np.int32)[0])
        L = float(np.frombuffer(f.read(8))[0]); Lc = float(np.frombuffer(f.read(8))[0])
        n = np.frombuffer(f.read(), np.float64).reshape(N, N, Nz, 3)
    hx, hz = L/(N-1), Lc/Nz
    d = [np.gradient(n, hx, axis=0), np.gradient(n, hx, axis=1), (np.roll(n, -1, 2) - np.roll(n, 1, 2))/(2*hz)]
    F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
    B = np.stack([F(1, 2), F(2, 0), F(0, 1)], -1)
    kx = 2*np.pi*np.fft.fftfreq(N, hx); kz = 2*np.pi*np.fft.fftfreq(Nz, hz)
    KX, KY, KZ = np.meshgrid(kx, kx, kz, indexing='ij'); K2 = KX**2 + KY**2 + KZ**2; K2[0, 0, 0] = 1
    Bk = [np.fft.fftn(B[..., i]) for i in range(3)]
    Ak = [1j*(KY*Bk[2] - KZ*Bk[1])/K2, 1j*(KZ*Bk[0] - KX*Bk[2])/K2, 1j*(KX*Bk[1] - KY*Bk[0])/K2]
    A = [np.real(np.fft.ifftn(x)) for x in Ak]
    Q = sum((A[i]*B[..., i]).sum() for i in range(3))*hx*hx*hz/(16*np.pi**2)
    flujo = B[..., 2].sum(axis=(0, 1)).mean()*hx*hx/(4*np.pi)   # grado de cada corte z (debe ser 0)
    print(f"{fich:30s} Lc={Lc:7.3f}  Q = {Q:+.4f}   grado de los cortes z = {flujo:+.4f}")
