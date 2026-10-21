# hopfq.py — carga de Hopf Q = (1/16π²) ∫ A·B d³x  (B_i = ½ ε_ijk F_jk, ∇×A = B, gauge de Coulomb por FFT)
import numpy as np, sys
f = open(sys.argv[1], 'rb'); N = np.frombuffer(f.read(4), np.int32)[0]; L = np.frombuffer(f.read(8))[0]
n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3); h = L/(N-1)
d = [np.gradient(n, h, axis=a) for a in range(3)]
F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
B = np.stack([F(1, 2), F(2, 0), F(0, 1)], -1)
k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij'); K2 = KX**2+KY**2+KZ**2; K2[0,0,0] = 1
Bk = [np.fft.fftn(B[..., i]) for i in range(3)]
# A = i k × B / k²
Ak = [1j*(KY*Bk[2]-KZ*Bk[1])/K2, 1j*(KZ*Bk[0]-KX*Bk[2])/K2, 1j*(KX*Bk[1]-KY*Bk[0])/K2]
A = [np.real(np.fft.ifftn(a)) for a in Ak]
Q = sum((A[i]*B[..., i]).sum() for i in range(3))*h**3/(16*np.pi**2)   # Q = (1/16π²)∫A·B
print(f"Q (carga de Hopf) = {Q:.4f}")
