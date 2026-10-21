// reposo2.cpp — Leptón en reposo en las dimensiones extendidas (simetría esférica, l = 0).
// Programa de referencia con el estado vigente (18-10-2026). († Claude, a partir de reposo.cpp del 11-10)
//
// QUÉ DESCRIBE
//   La envolvente en las extendidas, igual para electrón, muón y tau: misma ecuación, cada uno en sus
//   propias unidades (x en ħ/mc, energías en mc²). Lo que distingue a los tres leptones está en la
//   sección compacta (modos radiales de las capas A y C: columpio impulsado, sección 5 ter); aquí solo
//   cambia la conversión a fm. El muón y el tau ligados no tienen esta estructura: la envolvente se
//   convierte en el orbital.
//
// ECUACIONES
//   Núcleo (autopotencial V ∝ Ψ con estructura de Klein–Gordon; meseta en el umbral de Klein v = 2):
//       u'' + (2/x) u' = [ (1 + eps − a·u + w)² − 1 ] u ,   a = 1/[2(1+eps)] ,   eps = 0 en reposo
//   Carga: total Q = 1 (un cuanto), repartida según la opción peso=:
//       peso=0  ρ ∝ −∇²u = (1 − q²)u   (lectura vigente: Ψ como potencial, carga = −∇²Ψ;  [H])  <- por defecto
//               su integral es el flujo en el borde: Q ∝ −4π x1² u'(x1) (ley de Gauss)
//       peso=1  ρ ∝ |u'|               (agarre antiguo de reposo.cpp; superado, se deja para comparar)
//       peso=2  ρ ∝ u²·|u'|            (variante de prueba)
//       peso=3  ρ puntual (cuanto concentrado en las compactadas; con ell > 0 da el núcleo de Bopp–Podolsky)
//   Potencial eléctrico φ (Ψ como potencial; ley de Gauss), con regularización de Bopp–Podolsky
//   opcional (ell > 0):  (−∇² + ell²∇⁴) φ = 4πρ   ->   φ = φ_Coulomb[ρ] − φ_Yukawa,ell[ρ]
//       ell = 0: Coulomb puro.  <r²> efectivo = <r²>_ρ + 6 ell²  [D].  Cotas (muón): ⁴He -> ell ≲ 0,058 fm
//       (≈ 4α en unidades del muón); LEP -> r < 2,4e-4 fm, mucho más fuerte (pide el cierre dinámico) [A].
//   Acción del campo sobre el núcleo:  w = gw·φ.  gw = 0 POR DEFECTO (18-10): la envolvente es el estado del
//   leptón y su nube no puede repelerse a sí misma (el hidrógeno lo prohíbe) [D]. gw=α queda para comparar
//   (antes: energía de la carga en su propio
//   potencial, en unidades de mc²). Se itera hasta la autoconsistencia.
//
// USO
//   ./reposo2 u0_1 [u0_2 ...] [peso=0] [ell=0] [gw=0] [eps=0]
//   Proyectil externo (18-10, prueba de la compresión adiabática):  ext=s  b=b
//       w_ext = s/max(x,b): promedio esférico del potencial de una carga a distancia b (s = ±Z·α; + repele).
//       Limitación: con |s| grande el corte "vuelve a crecer" detiene la integración (x1 = 0,5 sin R50);
//       ese caso pide un tiro (shooting) en u0 que aún no está hecho.
//   PERFIL=1 ./reposo2 ...   escribe perfil2_u0_<u0>.dat :  x  u  rho  phi  w
// SALIDA (por u0): x1 (borde del núcleo, u = 0), R50 (radio a media altura), rq = <r²>^½ de la carga,
//   q_lejos (comprobación de Gauss: debe ser 1), E_campo = ½∫ρφ d³x en unidades de α·mc²,
//   y la conversión de R50 y rq a fm para electrón, muón y tau.
//
// ESTATUS (sinceridad)
//   [H] núcleo V ∝ Ψ con meseta de Klein (vigente).  [H] carga = −∇²Ψ.  [D] Gauss y Bopp–Podolsky.
//   [A] u0 es libre (no hay todavía condición que la fije); se barre.
//   [A] Con peso 0, 1 o 2 la carga ocupa toda la envolvente: rq ~ λ̄ (≈1600 fm para el electrón),
//       incompatible con la puntualidad medida de los tres leptones. Solo peso=3 (carga en las compactadas,
//       regla de las compactadas) es compatible; entonces la envolvente u es reparto de masa, no de carga.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
using namespace std;

const double ALFA = 1/137.035999;
const double LBAR[3] = {386.159, 1.86759, 0.111056};     // ħ/(mc) en fm: e, μ, τ
const char*  NOMBRE[3] = {"e", "mu", "tau"};

const double L = 60.0; const int N = 12000; const double dx = L/N;
vector<double> X(N);
vector<double> WEXT;   // potencial externo (proyectil), promedio esférico: s/max(x,b)

struct Nucleo { double x1 = NAN, R50 = NAN; vector<double> u, du; };

Nucleo integra(double u0, double eps, const vector<double>& w) {
    double a = 1.0/(2.0*(1.0+eps));
    auto F = [&](double x, double u, double du, double ww) {
        double q = 1.0 + eps - a*u + ww; return -2.0*du/x + (q*q - 1.0)*u; };
    Nucleo n; n.u.assign(N, 0.0); n.du.assign(N, 0.0);
    double q0 = 1.0 + eps - a*u0 + w[0] + WEXT[0], F0 = (q0*q0 - 1.0)*u0;
    double x = X[0], u = u0 + F0*x*x/6.0, du = F0*x/3.0;
    n.u[0] = u; n.du[0] = du;
    for (int i = 0; i < N-1; i++) {
        double w0 = w[i]+WEXT[i], w1 = w[i+1]+WEXT[i+1], wm = 0.5*(w0+w1);
        double k1u = du,              k1v = F(x, u, du, w0);
        double k2u = du + 0.5*dx*k1v, k2v = F(x+0.5*dx, u+0.5*dx*k1u, du+0.5*dx*k1v, wm);
        double k3u = du + 0.5*dx*k2v, k3v = F(x+0.5*dx, u+0.5*dx*k2u, du+0.5*dx*k2v, wm);
        double k4u = du + dx*k3v,     k4v = F(x+dx, u+dx*k3u, du+dx*k3v, w1);
        double un = u + dx/6.0*(k1u+2*k2u+2*k3u+k4u), dn = du + dx/6.0*(k1v+2*k2v+2*k3v+k4v);
        x += dx;
        if (std::isnan(n.R50) && un <= 0.5*u0) n.R50 = x;
        if (un <= 0.0) { n.x1 = x; break; }            // borde del núcleo
        if (dn > 0.0 && x > 0.5) { n.x1 = x; break; }   // vuelve a crecer: se corta aquí
        u = un; du = dn; n.u[i+1] = u; n.du[i+1] = du;
    }
    return n;
}

// densidad de carga normalizada a Q = 1 (∫ρ 4πx² dx = 1); peso 3 -> puntual (se trata aparte)
vector<double> densidad(const Nucleo& n, int peso, const vector<double>& w) {
    vector<double> rho(N, 0.0); double Q = 0;
    for (int i = 0; i < N; i++) {
        double q = 1.0 - 0.5*n.u[i] + w[i] + WEXT[i];           // eps = 0
        double wgt = (peso == 0) ? max(0.0, (1.0 - q*q)*n.u[i])
                   : (peso == 1) ? fabs(n.du[i]) : n.u[i]*n.u[i]*fabs(n.du[i]);
        rho[i] = wgt; Q += wgt*4*M_PI*X[i]*X[i]*dx;
    }
    for (double& r : rho) r = (Q > 0) ? r/Q : 0.0;
    return rho;
}

// potencial de Gauss (+ Bopp–Podolsky si ell > 0) de una densidad esférica
vector<double> potencial(const vector<double>& rho, int peso, double ell) {
    vector<double> phi(N, 0.0);
    if (peso == 3) {                                   // carga puntual en el origen
        for (int i = 0; i < N; i++) phi[i] = (ell > 0) ? (1 - exp(-X[i]/ell))/X[i] : 1/X[i];
        return phi;
    }
    // Coulomb: φ(x) = (1/x)∫0^x ρ 4πx'² dx' + ∫x^∞ ρ 4πx' dx'
    vector<double> qin(N), out(N);
    double s = 0; for (int i = 0; i < N; i++) { s += rho[i]*4*M_PI*X[i]*X[i]*dx; qin[i] = s; }
    s = 0; for (int i = N-1; i >= 0; i--) { out[i] = s; s += rho[i]*4*M_PI*X[i]*dx; }
    for (int i = 0; i < N; i++) phi[i] = qin[i]/X[i] + out[i];
    if (ell > 0) {                                     // resta del Yukawa (forma estable en exponenciales)
        double k = 1/ell;
        for (int i = 0; i < N; i++) {
            double r = X[i], sum = 0;
            for (int j = 0; j < N; j++) {
                double rp = X[j], rmin = min(r, rp), rmax = max(r, rp);
                double ker = 0.5*(exp(-k*(rmax-rmin)) - exp(-k*(rmax+rmin)));   // e^{-k r>} sinh(k r<)
                sum += rho[j]*rp*ker*dx;
            }
            phi[i] -= 4*M_PI*ell/r*sum;
        }
    }
    return phi;
}

int main(int argc, char** argv) {
    int peso = 0; double ell = 0, gw = 0, eps = 0, sext = 0, bext = 1; vector<double> u0s;
    for (int k = 1; k < argc; k++) {
        string a = argv[k];
        if (a.rfind("peso=", 0) == 0) peso = atoi(a.c_str()+5);
        else if (a.rfind("ell=", 0) == 0) ell = atof(a.c_str()+4);
        else if (a.rfind("ext=", 0) == 0) sext = atof(a.c_str()+4);
        else if (a.rfind("b=", 0) == 0) bext = atof(a.c_str()+2);
        else if (a.rfind("gw=", 0) == 0) gw = atof(a.c_str()+3);
        else if (a.rfind("eps=", 0) == 0) eps = atof(a.c_str()+4);
        else u0s.push_back(atof(a.c_str()));
    }
    if (u0s.empty()) { fprintf(stderr, "uso: %s u0_1 [u0_2 ...] [peso=0|1|2|3] [ell=0] [gw=0] [eps=0]\n", argv[0]); return 1; }
    if (peso == 3 && ell <= 0) fprintf(stderr, "aviso: carga puntual con ell=0 -> potencial 1/x singular en el centro\n");
    bool perfil = getenv("PERFIL") != nullptr;
    for (int i = 0; i < N; i++) X[i] = (i+1)*dx;
    WEXT.assign(N, 0.0);
    for (int i = 0; i < N; i++) WEXT[i] = sext/max(X[i], bext);

    printf("# reposo2: ext=%g b=%g peso=%d ell=%g gw=%g eps=%g   (x en ħ/mc, energías en mc²)\n", sext, bext, peso, ell, gw, eps);
    printf("#   u0      x1      R50      rq    q_lejos  E_campo[α·mc²] it |  R50[fm]: e  mu  tau  |  rq[fm]: e  mu  tau\n");
    for (double u0 : u0s) {
        vector<double> w(N, 0.0), rho, phi; Nucleo n; int it; double dw = 0;
        for (it = 0; it < 200; it++) {
            n = integra(u0, eps, w);
            rho = densidad(n, peso, w);
            phi = potencial(rho, peso, ell);
            dw = 0;
            for (int i = 0; i < N; i++) { double wn = gw*phi[i]; dw = max(dw, fabs(wn - w[i])); w[i] = 0.5*w[i] + 0.5*wn; }
            if (dw < 1e-9 || gw == 0) break;
        }
        double r2 = 0, E = 0;
        if (peso == 3) r2 = (ell > 0) ? 6*ell*ell : 0;   // Bopp–Podolsky: <r²> = 6 ell²
        else { for (int i = 0; i < N; i++) r2 += rho[i]*X[i]*X[i]*4*M_PI*X[i]*X[i]*dx; r2 += 6*ell*ell; }
        if (peso == 3) E = (ell > 0) ? 1/(2*ell) : INFINITY;
        else for (int i = 0; i < N; i++) E += 0.5*rho[i]*phi[i]*4*M_PI*X[i]*X[i]*dx;
        double rq = sqrt(r2), ql = phi[N-1]*X[N-1];
        printf("%8.4f %7.3f %8.4f %8.4f %8.4f %10.4f %3d |", u0, n.x1, n.R50, rq, ql, E, it);
        for (int l = 0; l < 3; l++) printf(" %9.4g", n.R50*LBAR[l]);
        printf("  |");
        for (int l = 0; l < 3; l++) printf(" %9.4g", rq*LBAR[l]);
        printf("\n"); fflush(stdout);
        if (perfil) {
            char nom[64]; snprintf(nom, sizeof nom, "perfil2_u0_%g.dat", u0);
            FILE* f = fopen(nom, "w");
            fprintf(f, "# x[ħ/mc]  u  rho  phi  w     (x en fm: multiplicar por %g (e), %g (mu), %g (tau))\n",
                    LBAR[0], LBAR[1], LBAR[2]);
            for (int i = 0; i < N; i += 10) fprintf(f, "%g %g %g %g %g\n", X[i], n.u[i], rho[i], phi[i], w[i]);
            fclose(f);
        }
    }
    printf("# Nota: los tres leptones son la misma solución; solo cambia la conversión a fm.\n");
    printf("# Ligados (átomos), muón y tau no tienen esta estructura: la envolvente pasa a ser el orbital.\n");
    return 0;
}
