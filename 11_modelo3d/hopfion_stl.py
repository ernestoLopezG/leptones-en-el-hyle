"""Modelo 3D (STL) de las preimágenes de un nudo de carga de Hopf 1: toros anidados como muñecas rusas,
cada uno cubierto por anillos de Villarceau enlazados, más el anillo central (el «alma» del nudo).
Geometría exacta: fibras de la fibración de Hopf proyectadas estereográficamente desde S³.

Uso: python3 hopfion_stl.py [capas=3] [anillos=8] [tamano_mm=100]
Salida: hopfion_todo.stl (todo junto) y hopfion_capa1.stl ... (cada capa por separado, para imprimir en colores)
Impresión: los anillos están enlazados pero no se tocan; hay que imprimirlos en polvo (SLS, nailon) o en resina
con soportes. En impresoras de filamento (FDM) conviene imprimir cada capa por separado y no se pueden enlazar luego.
"""
import sys, struct, numpy as np

capas  = int(sys.argv[1]) if len(sys.argv) > 1 else 3
anillos = int(sys.argv[2]) if len(sys.argv) > 2 else 8
tam_mm = float(sys.argv[3]) if len(sys.argv) > 3 else 100.0

def fibra(a, phi, n=400):
    """Fibra de Hopf: z1 = cos a·e^{it}, z2 = sin a·e^{i(t+phi)}, proyectada desde (0,0,0,1)."""
    t = np.linspace(0, 2*np.pi, n, endpoint=False)
    x1, x2 = np.cos(a)*np.cos(t), np.cos(a)*np.sin(t)
    x3, x4 = np.sin(a)*np.cos(t+phi), np.sin(a)*np.sin(t+phi)
    d = 1 - x4
    return np.stack([x1/d, x2/d, x3/d], axis=1)

def tubo(c, r, m=16):
    """Malla de un tubo cerrado de radio r alrededor de la curva cerrada c (transporte paralelo del marco)."""
    n = len(c); T = np.roll(c, -1, 0) - np.roll(c, 1, 0); T /= np.linalg.norm(T, axis=1)[:, None]
    v = np.cross(T[0], [0.3, 0.5, 0.8]); v /= np.linalg.norm(v); N = [v]
    for i in range(1, n):
        v = N[-1] - np.dot(N[-1], T[i])*T[i]; v /= np.linalg.norm(v); N.append(v)
    N = np.array(N)
    # corrige el giro acumulado para que el tubo cierre sin costura
    B0 = np.cross(T[-1], N[-1]); vn = N[0] - np.dot(N[0], T[-1])*T[-1]
    ang = np.arctan2(np.dot(vn, B0), np.dot(vn, N[-1]))
    s = np.arange(n)/n; cs, sn = np.cos(ang*s)[:, None], np.sin(ang*s)[:, None]
    B = np.cross(T, N); N, B = cs*N + sn*B, -sn*N + cs*B
    u = np.linspace(0, 2*np.pi, m, endpoint=False)
    P = c[:, None, :] + r*(np.cos(u)[None, :, None]*N[:, None, :] + np.sin(u)[None, :, None]*B[:, None, :])
    tri = []
    for i in range(n):
        for j in range(m):
            a, b = P[i, j], P[(i+1) % n, j]; cc, d = P[(i+1) % n, (j+1) % m], P[i, (j+1) % m]
            tri += [(a, b, cc), (a, cc, d)]
    return tri

def guarda(nombre, tri):
    with open(nombre, 'wb') as f:
        f.write(b'hopfion'.ljust(80, b' ')); f.write(struct.pack('<I', len(tri)))
        for a, b, c in tri:
            nrm = np.cross(b-a, c-a); l = np.linalg.norm(nrm); nrm = nrm/l if l > 0 else nrm
            f.write(struct.pack('<12fH', *nrm, *a, *b, *c, 0))
    print(f"{nombre}: {len(tri)} triángulos")

# capas: valores de a entre el alma (a=0) y el eje (a=pi/2)
As = [np.pi/2*(k+1)/(capas+2) for k in range(capas)]
curvas = {0: [fibra(0.0, 0.0)]}                      # el alma: circunferencia unidad
for k, a in enumerate(As, 1):
    curvas[k] = [fibra(a, 2*np.pi*j/anillos + 0.3*k) for j in range(anillos)]
todas = np.concatenate([c for cs in curvas.values() for c in cs])
escala = tam_mm/(todas.max(0) - todas.min(0)).max()
# radio del tubo: la mitad de la distancia mínima entre anillos distintos, con margen
pts = [c*escala for cs in curvas.values() for c in cs]
dmin = min(np.min(np.linalg.norm(pts[i][::4, None, :] - pts[j][None, ::4, :], axis=2))
           for i in range(len(pts)) for j in range(i+1, len(pts)))
r = min(0.35*dmin, 1.6)
print(f"tamaño {tam_mm:.0f} mm; distancia mínima entre anillos {dmin:.2f} mm; radio del tubo {r:.2f} mm")
total = []
for k, cs in curvas.items():
    tri = [t for c in cs for t in tubo(c*escala, r)]
    guarda(f"hopfion_capa{k}.stl" if k else "hopfion_alma.stl", tri); total += tri
guarda("hopfion_todo.stl", total)
np.save("curvas.npy", np.array([c*escala for cs in curvas.values() for c in cs]))
