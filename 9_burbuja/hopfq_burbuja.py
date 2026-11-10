# hopfq_burbuja.py — carga de Hopf del nudo dentro de la burbuja  (Claude, 10-11-2026)
# Q = (1/16π²) ∫ A·B d³x, con B_i = ½ ε_ijk F_jk y ∇×A = B (gauge de Coulomb por FFT), como hopfq.py,
# pero con B multiplicado por una ventana que apaga el fondo desordenado (fuera de él n no tiene sentido):
#   Q_completa : sin ventana (solo tiene sentido con fondo uniforme)
#   Q_S        : ventana = S(r)⁴, el mismo peso que el término de cuarto orden
#   Q_Rb       : ventana escalón en r < Rb + w
# Una ventana no da un entero exacto; lo que importa es si Q se queda cerca de ±1 o cae a 0.
# Uso: python3 hopfq_burbuja.py n_<etiqueta>.bin
import numpy as np, sys
f = open(sys.argv[1], 'rb')
N = int(np.frombuffer(f.read(4), np.int32)[0]); L, Rb, w = np.frombuffer(f.read(24))
n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3); h = L/(N-1)
d = [np.gradient(n, h, axis=a) for a in range(3)]
F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
B = np.stack([F(1, 2), F(2, 0), F(0, 1)], -1)
x = -L/2 + h*np.arange(N); X, Y, Z = np.meshgrid(x, x, x, indexing='ij'); r = np.sqrt(X**2 + Y**2 + Z**2)
S = 0.5*(1 - np.tanh((r - Rb)/w))
k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij'); K2 = KX**2 + KY**2 + KZ**2; K2[0, 0, 0] = 1
def Q(Bv):
    Bk = [np.fft.fftn(Bv[..., i]) for i in range(3)]
    Ak = [1j*(KY*Bk[2] - KZ*Bk[1])/K2, 1j*(KZ*Bk[0] - KX*Bk[2])/K2, 1j*(KX*Bk[1] - KY*Bk[0])/K2]
    A = [np.real(np.fft.ifftn(a)) for a in Ak]
    return sum((A[i]*Bv[..., i]).sum() for i in range(3))*h**3/(16*np.pi**2)
q_full = Q(B); q_S = Q(B*(S**4)[..., None]); q_R = Q(B*(r < Rb + w)[..., None])
print(f"{sys.argv[1]}  Rb={Rb:g}  Q_completa={q_full:.3f}  Q_S={q_S:.3f}  Q_Rb={q_R:.3f}")
