#!/usr/bin/env python3
# perturba.py — prueba de punto de silla para el hopfión (objeción M5 del árbitro)  (Claude, 7-11-2026)
#
# Lee un estado relajado de hopfion.cpp (n_final_N<N>.bin), lo deforma de forma NO simétrica y lo guarda
# para continuar la relajación con:  ./hopfion N L a iter perturbado.bin
#   deformación: estiramiento ×SX en x, ×SZ en z (rompe la simetría axial y la proporción del toro)
#                + ruido aleatorio de amplitud RUIDO en cada punto (rompe cualquier simetría restante)
# Si el nudo es un mínimo, al relajar vuelve a la misma energía y a la simetría axial.
# Si fuera un punto de silla, la energía bajaría a un valor distinto (y la forma no volvería).
#
# Uso: python3 perturba.py n_final_N81.bin  [SX=1.25] [SZ=0.8] [RUIDO=0.05]
# También: python3 perturba.py --asimetria fich.bin   → mide la asimetría de la forma
import sys, numpy as np
from scipy.ndimage import map_coordinates

def lee(fich):
    with open(fich, 'rb') as f:
        N = int(np.frombuffer(f.read(4), np.int32)[0]); L = float(np.frombuffer(f.read(8))[0])
        n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3).copy()
    return N, L, n

def asimetria(N, L, n):
    # momentos de segundo orden de w = 1 − n_z: Ixx, Iyy, Izz. Un toro de eje z tiene Ixx = Iyy.
    h = L/(N-1); x = -L/2 + h*np.arange(N)
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij'); w = 1 - n[..., 2]
    I = [(c**2*w).sum()/w.sum() for c in (X, Y, Z)]
    return I, abs(I[0] - I[1])/(I[0] + I[1])

if sys.argv[1] == '--asimetria':
    for fich in sys.argv[2:]:
        N, L, n = lee(fich); I, a = asimetria(N, L, n)
        print(f"{fich:40s} <x²>={I[0]:.4f} <y²>={I[1]:.4f} <z²>={I[2]:.4f}  asimetría x–y = {a:.4f}")
    sys.exit()

fich = sys.argv[1]
SX = float(sys.argv[2]) if len(sys.argv) > 2 else 1.25
SZ = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
RUIDO = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
N, L, n = lee(fich)
idx = np.arange(N); c = (N-1)/2
I, J, K = np.meshgrid(idx, idx, idx, indexing='ij')
coords = [c + (I - c)/SX, J, c + (K - c)/SZ]                  # n'(x,y,z) = n(x/SX, y, z/SZ)
m = np.stack([map_coordinates(n[..., p], coords, order=1, mode='nearest') for p in range(3)], -1)
rng = np.random.default_rng(1)
m += RUIDO*rng.standard_normal(m.shape)
m /= np.linalg.norm(m, axis=-1, keepdims=True)
m[0, :, :] = m[-1, :, :] = m[:, 0, :] = m[:, -1, :] = m[:, :, 0] = m[:, :, -1] = (0, 0, 1)
with open('perturbado.bin', 'wb') as f:
    np.array([N], np.int32).tofile(f); np.array([L]).tofile(f); m.astype(np.float64).tofile(f)
I0, a0 = asimetria(N, L, n); I1, a1 = asimetria(N, L, m)
print(f"original:   asimetría x–y = {a0:.4f}  <z²> = {I0[2]:.4f}")
print(f"perturbado: asimetría x–y = {a1:.4f}  <z²> = {I1[2]:.4f}   (SX={SX}, SZ={SZ}, ruido={RUIDO})")
print("-> perturbado.bin ; continuar con: ./hopfion", N, L, "1.0 <iter> perturbado.bin")
