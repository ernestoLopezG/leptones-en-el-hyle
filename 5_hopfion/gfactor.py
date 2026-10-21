# gfactor.py — factor g del hopfión con circulación compacta que avanza en el tiempo  († Claude, 20-10-2026)
#   Planteamiento [H]: n(x, η, t) = R_ẑ(kη − ωt)·n₃(x)  (la circulación compacta avanza a lo largo de η con el tiempo)
#   Pulsaciones sin masa en el hyle (Ernesto): avanzan a c a lo largo de η  ->  ω = k (c = 1).
#   Cuentas [D] con la energía de Faddeev–Skyrme en 4+1 (Lagrangiano ∂n·∂n − ½F·F):
#     ∂_t n = −ω ẑ×n , ∂_η n = k ẑ×n , F_ti = ω ∂_i n_z , F_iη = −k ∂_i n_z , F_tη = 0
#     I1 = ∫(1 − n_z²) ,  I2 = ∫|∇n_z|² ,  I = I1 + I2 ,  E_FS = ∫[∂n·∂n + ½F_ijF_ij]
#     Ecuación de n₃: extremo de E_FS + (k² − ω²)·I  ->  con ω = k es EXACTAMENTE el hopfión de Faddeev–Skyrme.
#     masa        M   = E_FS + (ω² + k²)·I = E_FS + 2ω²I        (T_00 = ṅ² + …, sin factor 2)
#     carga (KK)  Qe ∝ P_η = 2ωk·I            (impulso canónico a lo largo de η: T⁰_η = 2ṅ·∂_η n + …)
#     espín       J_z = 2m·ω·I                (T⁰_i = 2ṅ·∂_i n + …; m = vueltas alrededor de z)
#     mom. magn.  μ_z ∝ k·m·I                 (corriente de KK, T^{iη}, también con factor 2)
#     g = (2M/Qe)(μ/J) = M/(2ω²I) = 1 + E_FS/(2ω²·I)
#   CORRECCIÓN 20-10 9:30: la versión anterior omitía el factor 2 de las densidades de impulso y daba g = 2 (FALSO).
#   g = 1 es el de un momento orbital: un director (vector) que gira da espín entero. g = 2 y espín ½ exigen un espinor.
import numpy as np, sys
f = open(sys.argv[1], 'rb'); N = np.frombuffer(f.read(4), np.int32)[0]; L = np.frombuffer(f.read(8))[0]
n = np.frombuffer(f.read(), np.float64).reshape(N, N, N, 3); h = L/(N-1); dV = h**3
d = [np.gradient(n, h, axis=a) for a in range(3)]
F = lambda p, q: np.einsum('...k,...k->...', n, np.cross(d[p], d[q]))
E2 = sum((d[a]**2).sum() for a in range(3))*dV
E4 = sum((F(p, q)**2).sum() for p, q in [(0, 1), (0, 2), (1, 2)])*dV
I1 = (1 - n[..., 2]**2).sum()*dV
I2 = sum((d[a][..., 2]**2).sum() for a in range(3))*dV
EFS, I = E2 + E4, I1 + I2
print(f"E_FS = {EFS:.2f}  (E2 = {E2:.2f}, E4 = {E4:.2f})   I1 = {I1:.2f}   I2 = {I2:.2f}   I = {I:.2f}   E_FS/I = {EFS/I:.3f}")
for w in [1, 2, 5, 10, 30, 100]:
    print(f"  ω = k = {w:5g}  ->  g = 1 + {EFS/(2*w*w*I):.6f}")
print("Con un director que gira, g → 1 (orbital). El g = 2 exige que la circulación compacta sea espinorial (orden ½).")
