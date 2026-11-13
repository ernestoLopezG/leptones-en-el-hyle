"""¿Elipses u óvalos? Ajusta una elipse alineada con los ejes a cada capa n_z = const del nudo relajado, en el corte y = 0.
Uso: python3 forma_capas.py n_<tag>.bin [más archivos]
Columnas: capa n_z; radio del centro; semiejes radiales hacia dentro y hacia fuera; semieje según el eje;
asimetría (fuera − dentro)/(fuera + dentro) (0 en una elipse); residuo RMS del ajuste de elipse (% del radio);
b/a = semieje según el eje / semieje radial medio.
"""
import sys, struct, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def carga(f):
    b = open(f, 'rb').read(); N = struct.unpack('i', b[:4])[0]; L = struct.unpack('d', b[4:12])[0]
    return N, L, np.frombuffer(b[12:], dtype=np.float64).reshape(N, N, N, 3)
for f in sys.argv[1:]:
    N, L, n = carga(f); x = np.linspace(-L/2, L/2, N); j0 = N//2
    c = 0.5*(n[:, j0-1, :, 2] + n[:, j0, :, 2])
    print(f"{f}  (N = {N}, caja {L})")
    print("  n_z   centro   a_dentro  a_fuera   b_eje   asimetría  residuo%   b/a")
    for lev in [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6]:
        cs = plt.contour(x, x, c.T, levels=[lev]); plt.close()
        segs = [s for s in cs.allsegs[0] if len(s) > 20 and s[:, 0].mean() > 0]
        if not segs: continue
        P = max(segs, key=len); r, z = P[:, 0], P[:, 1]
        rtop = r[np.argmax(z)]; a_in, a_out = rtop - r.min(), r.max() - rtop; b = 0.5*(z.max() - z.min())
        M = np.stack([r*r, z*z, r, z], 1); A, B, C, D = np.linalg.lstsq(M, np.ones_like(r), rcond=None)[0]
        r0, z0 = -C/(2*A), -D/(2*B); K = 1 + A*r0*r0 + B*z0*z0; ae, be = np.sqrt(K/A), np.sqrt(K/B)
        res = np.sqrt(np.mean((np.hypot((r-r0)/ae, (z-z0)/be) - 1)**2))*100
        print(f"  {lev:5.1f}  {r0:6.3f}   {a_in:7.3f}  {a_out:7.3f}  {b:7.3f}   {(a_out-a_in)/(a_out+a_in):7.3f}   {res:6.2f}   {b/((a_in+a_out)/2):6.3f}")
