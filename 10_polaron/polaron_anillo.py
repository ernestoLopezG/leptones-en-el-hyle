"""Juguete: onda en un anillo (eta en [0,2pi)), antiperiódica, que genera su propio
patrón en un medio con inercia M (polarón). Sector de dos zonas: componentes k = 1/2 + 2n.
Ecuaciones (hbar = 1, radio 1, masa de la onda 1):
  onda:   i dPsi/dt = -1/2 d2Psi/deta2 - lam*phi*Psi
  medio:  M d2phi/dt2 = -kap*phi + D d2phi/deta2 + lam*|Psi|^2
(acoplamiento derivado de una sola energía, así que E y P se conservan)
Soluciones que derivan a velocidad angular v: Psi = e^{-iEt} u(eta - v t), phi = f(eta - v t).
Momento conservado (Noether): P = P_onda + P_medio,
  P_onda = <u|-i d/deta|u>,  P_medio = M v Int f'^2.
Uso: python3 polaron_anillo.py
"""
import numpy as np

N = 40                      # componentes k = 1/2 + 2n, n = -N..N
n = np.arange(-N, N+1); k = 0.5 + 2*n
J = 512; eta = np.linspace(0, 2*np.pi, J, endpoint=False)
E_k = np.exp(1j*np.outer(eta, k))/np.sqrt(2*np.pi)   # base en la malla

def resolver(lam, kap, M, D, v, c0=None, it=400):
    c = np.zeros(len(k), complex); c[N] = 1; c[N-1] = 0.3; c[N+1] = 0.3
    if c0 is not None: c = c0.copy()
    c /= np.linalg.norm(c)
    for _ in range(it):
        u = E_k @ c; rho = np.abs(u)**2
        rm = np.fft.rfft(rho)/J; m = np.arange(len(rm))
        den = kap + (D - M*v*v)*m**2
        f = np.fft.irfft(lam*rm/den*J, J)
        V = -lam*f
        Vm = np.fft.fft(V)/J                         # V(eta) = sum Vm[m] e^{i m eta}
        dk = (n[:, None] - n[None, :])*2            # diferencia de k (par)
        H = Vm[dk % J] + np.diag(0.5*(k - v)**2)
        w, vec = np.linalg.eigh(H)
        cn = vec[:, 0]
        cn *= np.exp(-1j*np.angle(np.vdot(c, cn)))
        c = 0.5*c + 0.5*cn; c /= np.linalg.norm(c)
    u = E_k @ c; rho = np.abs(u)**2
    rm = np.fft.rfft(rho)/J; m = np.arange(len(rm))
    f = np.fft.irfft(lam*rm/(kap + (D - M*v*v)*m**2)*J, J)
    fp = np.fft.irfft(1j*m*np.fft.rfft(f), J)
    Pw = float(np.sum(np.abs(c)**2*k))
    Pm = M*v*np.sum(fp**2)*2*np.pi/J
    contraste = (rho.max()-rho.min())/(rho.max()+rho.min())
    # energía en el laboratorio: onda + medio (cinética y potencial) + acoplamiento
    Ew = np.sum(np.abs(c)**2*0.5*k**2)
    Ec = -lam*np.sum(f*rho)*2*np.pi/J
    Emed = (0.5*kap*np.sum(f**2) + 0.5*D*np.sum(fp**2) + 0.5*M*v*v*np.sum(fp**2))*2*np.pi/J
    return dict(c=c, Pw=Pw, Pm=Pm, P=Pw+Pm, contraste=contraste, E=Ew+Ec+Emed)

def v_para_medio(lam, kap, M, D, objetivo=0.5):
    """busca la deriva v con P_total = objetivo (bisección)"""
    s0 = resolver(lam, kap, M, D, 0.0)
    vmax = np.sqrt(D/M)*0.999 if M > 0 else 50.0   # deriva subsónica en el medio
    if resolver(lam, kap, M, D, vmax, s0['c'], it=150)['P'] < objetivo: return None, s0, None
    lo, hi = 0.0, vmax
    for _ in range(30):
        mid = 0.5*(lo+hi)
        if resolver(lam, kap, M, D, mid, s0['c'], it=150)['P'] < objetivo: lo = mid
        else: hi = mid
    return hi, s0, resolver(lam, kap, M, D, hi, s0['c'])

if __name__ == '__main__':
    kap, D = 1.0, 0.2
    print("lam  M      contraste  P_onda(v=0)  v_deriva    P_onda   P_medio   P_total  dE(deriva)")
    for lam in [np.sqrt(10.0), np.sqrt(12.0), np.sqrt(16.0)]:
        for M in [1.0, 10.0, 100.0, 1000.0, 10000.0]:
            v, s0, s = v_para_medio(lam, kap, M, D)
            if s is None:
                print(f"{lam*lam:4.0f} {M:6.0f}  sin solución"); continue
            print(f"{lam*lam:4.0f} {M:6.0f}  {s0['contraste']:9.3f}  {s0['Pw']:10.4f}  {v:10.5f}  {s['Pw']:7.4f}  {s['Pm']:8.4f}  {s['P']:7.4f}  {s['E']-s0['E']:10.2e}")
