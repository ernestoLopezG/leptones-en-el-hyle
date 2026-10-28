#!/usr/bin/env python3
# mathieu_q.py — funciones angulares de Mathieu para q negativo grande  († Claude, 23-10-2026)
# Ecuación: y'' + (a − 2q·cos 2η) y = 0.  Con q < 0 el potencial 2q·cos 2η tiene sus pozos en η = 0 y η = π.
# Solución de Floquet de orden ν: y = Σ_k A_k · e^{i(2k+ν)η}  →  (a − (2k+ν)²) A_k − q (A_{k−1} + A_{k+1}) = 0  (método de Hill)
# ν = 0 → ce₀ (periódica, periodo π);  ν = ½ → la del modelo (antiperiódica en una vuelta: y(η+2π) = −y(η)).
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def floquet(q, nu, K=None):
    K = K or int(max(40, 3*np.sqrt(abs(q)) + 40))
    k = np.arange(-K, K+1); M = np.diag((2*k + nu)**2).astype(float)
    M += q*(np.eye(len(k), k=1) + np.eye(len(k), k=-1))      # a·A = M·A
    w, v = np.linalg.eigh(M); i = 0                            # el valor propio más bajo (estado fundamental)
    return w[i], k, v[:, i]
eta = np.linspace(0, 2*np.pi, 2001)
qs = [-10, -100, -1000]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
for q in qs:
    a, k, A = floquet(q, 0.0)
    y = np.real(np.exp(1j*np.outer(eta, 2*k)) @ A); y /= np.max(np.abs(y))
    ax[0].plot(np.degrees(eta), y, label=f'q = {q}  (a = {a:.1f})')
    a, k, A = floquet(q, 0.5)
    psi = np.exp(1j*np.outer(eta, 2*k + 0.5)) @ A               # onda que circula en un solo sentido
    re = psi.real.copy(); re /= np.max(np.abs(re))   # copia: np.real devuelve una vista
    ax[1].plot(np.degrees(eta), re, label=f'q = {q}')
    I = np.abs(psi)**2; I /= I.max()
    ax[2].plot(np.degrees(eta), I, label=f'q = {q}')
for a_ in ax:
    a_.set_xlabel('η (grados)'); a_.set_xticks(range(0, 361, 90)); a_.grid(alpha=.3); a_.legend(fontsize=8)
    a_.axvline(0, c='k', lw=.5); a_.axvline(180, c='k', lw=.5)
ax[0].set_title('ce₀(η): orden entero (periodo π)')
ax[1].set_title('Orden ν = ½: parte real de la onda que circula\n(cambia de signo tras una vuelta)')
ax[2].set_title('Intensidad |Ψ|² de la onda de orden ½\n(dos máximos opuestos: las nucleaciones)')
fig.suptitle('Funciones angulares de Mathieu con q negativo grande: se concentran en los pozos η = 0° y η = 180°', fontsize=11)
plt.tight_layout(); plt.savefig('mathieu_q_negativo.png', dpi=120); print('-> mathieu_q_negativo.png')
for q in qs:
    a, k, A = floquet(q, 0.5); a0, *_ = floquet(q, 0.0)
    print(f'q = {q:6d}: a(ν=½) = {a:.6f}   a(ν=0) = {a0:.6f}   diferencia = {a-a0:.3e}')
