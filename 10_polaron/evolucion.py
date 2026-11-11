"""Evolución temporal del juguete del polarón en el anillo, para comprobar que:
 (i) el momento total P = P_onda + P_medio se conserva;
 (ii) el estado que deriva (calculado con polaron_anillo.py) sigue derivando a la
      velocidad v prevista, con el patrón de dos zonas intacto.
Onda: paso partido (cinética en Fourier, potencial en la malla). Medio: Verlet.
Uso: python3 evolucion.py lam M  (por defecto lam=sqrt(10), M=100)
"""
import sys, numpy as np, importlib.util
spec = importlib.util.spec_from_file_location('p', 'polaron_anillo.py')
p = importlib.util.module_from_spec(spec); spec.loader.exec_module(p)

lam = float(sys.argv[1]) if len(sys.argv) > 1 else np.sqrt(10.0)
M = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
kap, D = 1.0, 0.2
v, s0, s = p.v_para_medio(lam, kap, M, D)
print(f"lam^2={lam*lam:.1f} M={M:g}: v prevista={v:.5f}, P_onda={s['Pw']:.4f}, P_medio={s['Pm']:.4f}")

J = p.J; eta = p.eta; dx = 2*np.pi/J
u = p.E_k @ s['c']                                   # onda inicial
rho = np.abs(u)**2
mm = np.fft.fftfreq(J, 1/J)
f = np.real(np.fft.ifft(lam*np.fft.fft(rho)/(kap + (D - M*v*v)*mm**2)))
fp = np.real(np.fft.ifft(1j*mm*np.fft.fft(f)))
phi, pi_ = f.copy(), M*v*fp                          # phi_t = -v f'  ->  momento M phi_t; aquí pi_ = -M phi_t
phit = -v*fp
# la onda es antiperiódica: se trabaja con w = e^{-i eta/2} u, periódica
w = np.exp(-0.5j*eta)*u
kk = mm + 0.5                                         # número de onda real de cada componente

def momentos(w, phi, phit):
    wk = np.fft.fft(w)/J
    Pw = np.sum(np.abs(wk)**2*kk)/np.sum(np.abs(wk)**2)
    phix = np.real(np.fft.ifft(1j*mm*np.fft.fft(phi)))
    Pm = -M*np.sum(phit*phix)*dx
    return Pw, Pm

def lap(a): return np.real(np.fft.ifft(-(mm**2)*np.fft.fft(a)))

def centro(rho):                                      # ángulo del patrón de periodo pi
    z = np.sum(rho*np.exp(2j*eta)); return np.angle(z)/2

dt = 2e-4; T = 20.0 if M >= 100 else 4.0
nsteps = int(T/dt); salida = max(1, nsteps//10)
fase_k = np.exp(-0.5j*kk**2*dt)
a0 = centro(np.abs(w)**2); t = 0.0; ang_prev = a0; ang_acum = 0.0
print("   t      P_onda    P_medio   P_total   giro del patrón   contraste")
for i in range(nsteps+1):
    if i % salida == 0:
        Pw, Pm = momentos(w, phi, phit); r = np.abs(w)**2
        c = (r.max()-r.min())/(r.max()+r.min())
        print(f"{t:7.3f}  {Pw:8.4f}  {Pm:8.4f}  {Pw+Pm:8.4f}  {ang_acum:12.5f}  (v*t={v*t:.5f})  {c:6.3f}")
    # medio: Verlet
    acc = (-kap*phi + D*lap(phi) + lam*np.abs(w)**2)/M
    phit_half = phit + 0.5*dt*acc
    phi = phi + dt*phit_half
    # onda: paso partido de Strang con el potencial en el punto medio
    V = -lam*phi
    w = w*np.exp(-0.5j*V*dt)
    w = np.fft.ifft(fase_k*np.fft.fft(w))
    w = w*np.exp(-0.5j*V*dt)
    acc = (-kap*phi + D*lap(phi) + lam*np.abs(w)**2)/M
    phit = phit_half + 0.5*dt*acc
    t += dt
    a = centro(np.abs(w)**2); d = (a - ang_prev + np.pi/2) % np.pi - np.pi/2
    ang_acum += d; ang_prev = a
