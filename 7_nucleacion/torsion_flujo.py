# torsion_flujo.py — torsión del flujo inducido de nucleacion.cpp  (Claude, 2-11-2026)
#
# El flujo es un anillo de vórtice con remolino: Ω = g(d)·(φ̂ + s·ê_pol), g = exp(−d²/w²),
# d = distancia al núcleo del anillo (radio R0). Su helicidad es H = ∫ u·Ω con ∇×u = Ω.
# Para un tubo de vórtice, H = T·Φ², con Φ el flujo a lo largo del tubo y T su torsión: el número de
# vueltas que dan las líneas alrededor del núcleo por cada vuelta alrededor del eje (Moffatt y Ricca 1992).
# Hipótesis a contrastar con el barrido: la carga de Hopf del nudo atado es Q ≈ T redondeado.
#
# Uso: python3 torsion_flujo.py            (imprime T para varios s y guarda torsion_flujo.png)
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def campo(s, N=128, L=8.0, R0=1.2, w=0.6):
    x = np.linspace(-L/2, L/2, N, endpoint=False); h = x[1]-x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho = np.sqrt(X**2+Y**2) + 1e-12
    d = np.sqrt((rho-R0)**2 + Z**2) + 1e-12; g = np.exp(-d**2/w**2)
    er = np.stack([X/rho, Y/rho, 0*X]); ephi = np.stack([-Y/rho, X/rho, 0*X]); ez = np.stack([0*X, 0*X, 1+0*X])
    epol = (-Z*er + (rho-R0)*ez)/d
    return g*(ephi + s*epol), h, N

def helicidad(Om, h, N):
    k = 2*np.pi*np.fft.fftfreq(N, h); KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX**2+KY**2+KZ**2; K2[0, 0, 0] = 1
    Ok = [np.fft.fftn(c) for c in Om]
    div = (KX*Ok[0]+KY*Ok[1]+KZ*Ok[2])/K2              # parte solenoidal de Ω (quita la pequeña divergencia)
    Ok = [Ok[0]-KX*div, Ok[1]-KY*div, Ok[2]-KZ*div]
    Uk = [1j*(KY*Ok[2]-KZ*Ok[1])/K2, 1j*(KZ*Ok[0]-KX*Ok[2])/K2, 1j*(KX*Ok[1]-KY*Ok[0])/K2]   # u = ∇×Ω/(−∇²)... ∇×u = Ω
    u = [np.real(np.fft.ifftn(a)) for a in Uk]; O = [np.real(np.fft.ifftn(a)) for a in Ok]
    return sum((u[i]*O[i]).sum() for i in range(3))*h**3

if __name__ == "__main__":
    w = 0.6; Phi = np.pi*w**2                            # flujo a lo largo del tubo: ∫∫ g dA = π w²
    ss = np.array([0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0])
    T = []
    for s in ss:
        Om, h, N = campo(s); H = helicidad(Om, h, N); T.append(H/Phi**2)
        print(f"s = {s:5.2f}   H = {H:8.4f}   torsión T = H/Φ² = {H/Phi**2:6.3f}")
    T = np.array(T)
    obs = {0: [0, 0, 0], 0.3: [1]*7+[2], 0.5: [2]*6+[1], 1.0: [3, 3, 3]}   # barridos 1 y 2 (|Q| final)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ss, np.abs(T), "k-", label="torsión del flujo |T| = |H|/Φ²")
    for s, qs in obs.items():
        vals, cnt = np.unique(qs, return_counts=True)
        for v, c in zip(vals, cnt): ax.scatter([s], [v], s=25*c, color="C3", alpha=0.7)
    ax.scatter([], [], color="C3", label="|Q| del nudo atado (tamaño ∝ n.º de semillas)")
    for q in (1, 2, 3): ax.axhline(q, color="0.85", lw=0.6)
    ax.set_xlabel("remolino |s|"); ax.set_ylabel("vueltas / carga"); ax.legend(fontsize=8)
    ax.set_title("¿Hereda el nudo la torsión del flujo de la pulsación?", fontsize=10)
    fig.tight_layout(); fig.savefig("torsion_flujo.png", dpi=120); print("guardada torsion_flujo.png")
