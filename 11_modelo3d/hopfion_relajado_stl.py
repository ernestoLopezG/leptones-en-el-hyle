"""STL de las preimágenes del nudo RELAJADO (simulado), para compararlo con el mapa de Hopf ideal (hopfion_stl.py).
Lee un campo del director guardado por 7_nucleacion/nucleacion_pulsos (n_<tag>.bin: int N, double L, N³ vectores),
extrae los anillos donde el director apunta a una dirección dada (capas n_z = -0,6, 0 y 0,6 y el alma) y los suaviza
con pocos armónicos del ángulo alrededor del eje. Escala: el alma mide lo mismo que en el modelo ideal.

Uso: python3 hopfion_relajado_stl.py n_base62b.bin [anillos=6] [mm_por_unidad=auto]
"""
import sys, struct, numpy as np
from hopfion_stl import tubo, guarda   # reutiliza la malla de tubos (no ejecuta el modelo ideal: ver más abajo)

f = sys.argv[1]; anillos = int(sys.argv[2]) if len(sys.argv) > 2 else 6
b = open(f, 'rb').read(); N = struct.unpack('i', b[:4])[0]; L = struct.unpack('d', b[4:12])[0]
n = np.frombuffer(b[12:], dtype=np.float64).reshape(N, N, N, 3)
x = np.linspace(-L/2, L/2, N); X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

def anillo(c, phi, tol=0.15, nb=48, armonicos=3):
    s = np.sqrt(max(0.0, 1 - c*c)); t = np.array([s*np.cos(phi), s*np.sin(phi), c])
    m = np.linalg.norm(n - t, axis=-1) < tol
    P = np.stack([X[m], Y[m], Z[m]], 1); w = 1 - np.linalg.norm(n[m] - t, axis=-1)/tol
    psi = np.arctan2(P[:, 1], P[:, 0]); rho = np.hypot(P[:, 0], P[:, 1])
    # ajuste por mínimos cuadrados de rho(psi) y z(psi) con pocos armónicos
    cols = [np.ones_like(psi)] + [f(k*psi) for k in range(1, armonicos+1) for f in (np.cos, np.sin)]
    A = np.stack(cols, 1)*w[:, None]
    cr = np.linalg.lstsq(A, rho*w, rcond=None)[0]; cz = np.linalg.lstsq(A, P[:, 2]*w, rcond=None)[0]
    ps = np.linspace(-np.pi, np.pi, 400, endpoint=False)
    B = np.stack([np.ones_like(ps)] + [f(k*ps) for k in range(1, armonicos+1) for f in (np.cos, np.sin)], 1)
    r_, z_ = B@cr, B@cz
    return np.stack([r_*np.cos(ps), r_*np.sin(ps), z_], 1)

capas = [-0.6, 0.0, 0.6]
curvas = {0: [anillo(-0.97, 0.0, tol=0.2)]}
for k, c in enumerate(capas, 1):
    curvas[k] = [anillo(c, 2*np.pi*j/anillos + 0.3*k) for j in range(anillos)]
Rcore = np.hypot(curvas[0][0][:, 0], curvas[0][0][:, 1]).mean()
mm = float(sys.argv[3]) if len(sys.argv) > 3 else 120/ (2*np.abs(np.concatenate([c for cs in curvas.values() for c in cs])).max())
print(f"alma: radio {Rcore:.3f} unidades; escala {mm:.1f} mm por unidad")
pts = [c*mm for cs in curvas.values() for c in cs]
dmin = min(np.min(np.linalg.norm(pts[i][::4, None, :] - pts[j][None, ::4, :], axis=2))
           for i in range(len(pts)) for j in range(i+1, len(pts)))
r = min(0.35*dmin, 1.6); print(f"distancia mínima entre anillos {dmin:.2f} mm; radio del tubo {r:.2f} mm")
total = []
for k, cs in curvas.items():
    tri = [t for c in cs for t in tubo(c*mm, r)]
    guarda(f"relajado_capa{k}.stl" if k else "relajado_alma.stl", tri); total += tri
guarda("relajado_todo.stl", total)
np.save("curvas_relajado.npy", np.array(pts))
