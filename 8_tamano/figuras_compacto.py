#!/usr/bin/env python3
# figuras_compacto.py — resultados de hopfion_compacto.cpp  (Claude, 4-11-2026)
# Uso:  python3 figuras_compacto.py compacto_<etiqueta>.dat
# Dibuja: (a) tamaño en el plano extendido R_xy y anchura en z frente al periodo Lc;
#         (b) energía E/(32√2π²) frente a Lc;  (c) cortes de n_z del último y del primer estado (plano x–z por el centro).
import sys, glob, re, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
dat = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob('compacto_*.dat'))[0]
tag = re.sub(r'^compacto_|\.dat$', '', dat.split('/')[-1])
D = np.loadtxt(dat); D = D[np.argsort(-D[:, 0])]
Lc, E, Rxy, sz, Rbar, mm = D[:, 0], D[:, 4], D[:, 5], D[:, 6], D[:, 7], D[:, 8]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
ax[0].plot(Lc, Rxy, 'o-', label='R_xy (extendidas)')
ax[0].plot(Lc, np.where(np.isfinite(sz), sz*np.sqrt(2), np.nan), 's-', label='√2·σ_z (compacta)')
ax[0].plot(Lc, Lc/np.pi, 'k:', label='Lc/π (escala de la dirección compacta)')
ax[0].set_xlabel('periodo de la dirección compacta Lc'); ax[0].set_ylabel('tamaño')
ax[0].invert_xaxis(); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
ax[0].set_title('¿Sigue el tamaño extendido al compacto?')
a2 = ax[0].twinx(); a2.plot(Lc, mm, 'g^--', alpha=.6); a2.set_ylabel('min/max a lo largo de z (1 = llena la compacta)', color='g'); a2.set_ylim(0, 1.05)
ax[1].plot(Lc, E, 'o-'); ax[1].set_xlabel('Lc'); ax[1].set_ylabel('E/(32√2π²)'); ax[1].invert_xaxis(); ax[1].grid(alpha=.3)
ax[1].axhline(1.22, color='gray', ls=':', label='hopfión libre en R³ (≈1,22)'); ax[1].legend(fontsize=8); ax[1].set_title('Energía')
fichs = sorted(glob.glob(f'n_{tag}_Lc*.bin'), key=lambda s: -float(re.search(r'Lc([\d.]+)\.bin', s).group(1)))
if fichs:
    def lee(fich):
        with open(fich, 'rb') as f:
            N = int(np.frombuffer(f.read(4), np.int32)[0]); Nz = int(np.frombuffer(f.read(4), np.int32)[0])
            L = float(np.frombuffer(f.read(8))[0]); Lc_ = float(np.frombuffer(f.read(8))[0])
            n = np.frombuffer(f.read(), np.float64).reshape(N, N, Nz, 3)
        return N, Nz, L, Lc_, n
    N, Nz, L, Lc_, n = lee(fichs[-1])
    corte = n[:, N//2, :, 2]
    ax[2].imshow(corte.T, origin='lower', extent=[-L/2, L/2, 0, Lc_], aspect='equal', cmap='coolwarm', vmin=-1, vmax=1)
    ax[2].set_xlabel('x (extendida)'); ax[2].set_ylabel('z (compacta, periódica)')
    ax[2].set_title(f'n_z en el plano x–z, Lc = {Lc_:g}')
plt.tight_layout(); out = f'compacto_{tag}.png'; plt.savefig(out, dpi=120); print('->', out)
