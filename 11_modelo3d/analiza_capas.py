import numpy as np, struct
def carga(f):
    b=open(f,'rb').read(); N=struct.unpack('i',b[:4])[0]; L=struct.unpack('d',b[4:12])[0]
    n=np.frombuffer(b[12:],dtype=np.float64).reshape(N,N,N,3); return N,L,n
def hopf(N,L,a=1.0):
    x=np.linspace(-L/2,L/2,N); X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); X,Y,Z=X/a,Y/a,Z/a; r2=X*X+Y*Y+Z*Z
    u1,u2,u3,u4=2*X/(1+r2),2*Y/(1+r2),2*Z/(1+r2),(r2-1)/(1+r2)
    z1r,z1i,z0r,z0i=u1,u2,u4,u3; re=z1r*z0r+z1i*z0i; im=z1i*z0r-z1r*z0i; m1=z1r**2+z1i**2; m0=z0r**2+z0i**2
    return np.stack([2*re/(m0+m1),2*im/(m0+m1),(m0-m1)/(m0+m1)],-1)
N,L,nr=carga(__import__('sys').argv[1]); ni=hopf(N,L)
x=np.linspace(-L/2,L/2,N); X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); RHO=np.hypot(X,Y)
def toro(n,c,tol=0.06):
    m=np.abs(n[...,2]-c)<tol; rho,z=RHO[m],Z[m]
    # sección en el plano (rho,z): centro y semiejes
    R=rho.mean(); zc=z.mean(); a_r=np.sqrt(2)*np.std(rho); a_z=np.sqrt(2)*np.std(z)
    return R,zc,a_r,a_z,m.sum()
def preimagen(n,c,phi,tol=0.12,nb=36):
    t=np.array([np.sqrt(1-c*c)*np.cos(phi),np.sqrt(1-c*c)*np.sin(phi),c])
    m=np.linalg.norm(n-t,axis=-1)<tol; P=np.stack([X[m],Y[m],Z[m]],1)
    ang=np.arctan2(P[:,1],P[:,0]); bins=np.linspace(-np.pi,np.pi,nb+1); C=[]
    for i in range(nb):
        s=(ang>=bins[i])&(ang<bins[i+1])
        if s.sum(): C.append(P[s].mean(0))
    C=np.array(C); c0=C.mean(0); u,sv,vt=np.linalg.svd(C-c0); nrm=vt[2]
    plan=np.abs((C-c0)@nrm).max(); q=(C-c0)-np.outer((C-c0)@nrm,nrm); rad=np.linalg.norm(q,axis=1)
    return C, plan, rad.mean(), rad.std()/rad.mean(), np.degrees(np.arccos(abs(nrm[2])))
print("N",N,"L",L)
print("capa n_z     |  ideal: R    z0   semieje_rho semieje_z | relajado: R    z0   semieje_rho semieje_z")
for c in [-0.9,-0.6,-0.3,0.0,0.3,0.6]:
    a=toro(ni,c); b=toro(nr,c)
    print(f"{c:5.1f}  | {a[0]:6.3f} {a[1]:6.3f} {a[2]:7.3f} {a[3]:7.3f}   | {b[0]:6.3f} {b[1]:6.3f} {b[2]:7.3f} {b[3]:7.3f}")
print("\npreimágenes (anillos): planitud máx, radio, desviación relativa del radio, inclinación del plano")
for c in [-0.6,0.0,0.6]:
    for lab,n in [('ideal',ni),('relaj',nr)]:
        C,pl,rad,dev,inc=preimagen(n,c,0.0)
        print(f"n_z={c:4.1f} {lab}: planitud {pl:.3f}  radio {rad:.3f}  desv {dev:.3f}  inclinación {inc:.1f}°  ({len(C)} puntos)")

