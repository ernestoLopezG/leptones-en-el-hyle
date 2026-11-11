"""Solapamiento |<Psi(eta)|Psi(eta - pi/2)>| de dos ondas de Mathieu de orden 1/2 con las
líneas de zonas de mayor índice perpendiculares (espín opuesto en la imagen del preprint).
Método de Hill. Uso: python3 solapamientos.py"""
import numpy as np
def onda(q, N=60):
    n = np.arange(-N, N+1); k = 0.5 + 2*n
    M = np.diag(k**2) + q*(np.eye(len(k), k=1) + np.eye(len(k), k=-1))
    w, v = np.linalg.eigh(M); return k, v[:, 0]
for q in [0, -0.64, -1, -2, -3, -10, -100, -1000]:
    k, c = onda(q, N=60 if abs(q) < 200 else 200)
    print(f"q = {q:8}: solapamiento = {abs(np.sum(np.abs(c)**2*np.exp(-1j*k*np.pi/2))):.3g}")
