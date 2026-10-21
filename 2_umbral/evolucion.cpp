// evolucion.cpp — Estabilidad del solitón con umbral inferior (Bingham), evolución temporal radial († Claude, 18-10-2026)
//
// Ecuación de onda cuya forma estacionaria ψ = u(r)·e^{−it} es la de umbral.cpp:
//     ∂²ψ/∂t² − ∇²ψ + M²ψ = 0 ,   M = 1 − a·(|ψ| − uc)₊ ,   a = 1/2      (unidades: ħ = c = m = 1)
// Es Klein–Gordon con una "masa" que el propio hyle rebaja donde fluye (|ψ| > uc). En reposo, ω = 1:
//     ∇²u = (M² − 1)u   — exactamente la ecuación de umbral.cpp.  Fuera del núcleo M = 1: Klein–Gordon libre.
// [H] Esta forma temporal es la más sencilla compatible con la estática; no está derivada.
//
// Método: φ = r·ψ (complejo), φ_tt = φ_rr − M²φ, φ(0) = 0; diferencias centradas (salto de rana);
//   capa absorbente en el último tramo para que la radiación salga.
// Uso:  ./evolucion uc delta [tmax]
//   delta: perturbación inicial de la amplitud, ψ(0) = (1+delta)·u·e^{−it}  (delta = 0: prueba de estacionariedad)
// Salida: evol_uc<uc>_d<delta>.dat :  t  |ψ(0)|/u0  N(r<Rn)/N0  E... ; y en pantalla un resumen
//   N = ∫ Im(ψ* ψ_t) d³x (la "carga" conservada de la simetría de fase) dentro de r < 4R
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <complex>
using namespace std;
typedef complex<double> C;
const double a = 0.5;

// perfil estacionario por tiro (como umbral.cpp)
double uc_g;
double tiro(double u0, double& R, vector<double>* X = 0, vector<double>* U = 0) {
    double uc = uc_g, h = 1e-4;
    auto F = [&](double x, double u, double v) { double q = 1 - a*fmax(0.0, u - uc); return -2*v/x + (q*q - 1)*u; };
    double q0 = 1 - a*(u0 - uc), F0 = (q0*q0 - 1)*u0;
    double x = 1e-4, u = u0 + F0*x*x/6, v = F0*x/3;
    while (x < 500) {
        if (X) { X->push_back(x); U->push_back(u); }
        double k1u = v, k1v = F(x, u, v);
        double k2u = v + .5*h*k1v, k2v = F(x + .5*h, u + .5*h*k1u, v + .5*h*k1v);
        double k3u = v + .5*h*k2v, k3v = F(x + .5*h, u + .5*h*k2u, v + .5*h*k2v);
        double k4u = v + h*k3v,    k4v = F(x + h, u + h*k3u, v + h*k3v);
        double un = u + h/6*(k1u + 2*k2u + 2*k3u + k4u), vn = v + h/6*(k1v + 2*k2v + 2*k3v + k4v);
        if (un <= uc) { double t = (u - uc)/(u - un); R = x + t*h; return R*(v + t*(vn - v)) + uc; }
        if (vn > 0) { R = NAN; return NAN; }
        u = un; v = vn; x += h;
    }
    R = NAN; return NAN;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "uso: %s uc delta [tmax]\n", argv[0]); return 1; }
    uc_g = atof(argv[1]); double delta = atof(argv[2]), tmax = argc > 3 ? atof(argv[3]) : 600;
    double uc = uc_g;
    // u0 por bisección
    double prevm = NAN, prevu = NAN, u0 = NAN;
    for (int i = 1; i <= 4000 && std::isnan(u0); i++) {
        double uu = uc + 4.0*i/4001, R; double m = tiro(uu, R);
        if (!std::isnan(m) && !std::isnan(prevm) && m*prevm < 0) {
            double lo = prevu, hi = uu, mlo = prevm;
            for (int it = 0; it < 60; it++) { double mid = .5*(lo+hi), Rm; double mm = tiro(mid, Rm);
                if (std::isnan(mm)) break; if (mm*mlo < 0) hi = mid; else { lo = mid; mlo = mm; } }
            u0 = .5*(lo+hi);
        }
        prevm = m; prevu = uu;
    }
    double R; vector<double> X, U; tiro(u0, R, &X, &U);
    double c = uc*R;
    printf("uc = %g  u0 = %.6f  R = %.4f  c = %.4f  delta = %g\n", uc, u0, R, c, delta);

    // malla
    const double dr = getenv("DR") ? atof(getenv("DR")) : 0.02, L = 400, Labs = 120, dt = 0.4*dr;
    int N = (int)(L/dr) + 1;
    auto perfil = [&](double r) {               // u(r) interpolado + cola c/r
        if (r >= R) return c/r;
        double h = 1e-4; size_t j = (size_t)((r - 1e-4)/h); if (j + 1 >= X.size()) j = X.size() - 2;
        double t = (r - X[j])/h; return U[j] + t*(U[j+1] - U[j]); };
    vector<C> f(N), fo(N), fn(N);
    vector<double> sig(N, 0.0);
    for (int i = 0; i < N; i++) { double r = i*dr; if (r > L - Labs) { double s = (r - (L - Labs))/Labs; sig[i] = 2.0*s*s; } }
    // condición inicial: ψ = (1+δ) u e^{−it} ; paso atrás exacto de la fase: ψ(−dt) = (1+δ) u e^{+i dt}
    for (int i = 0; i < N; i++) { double r = i*dr, u = (i ? perfil(r) : 0.0)*(1 + delta);
        f[i] = r*u; fo[i] = r*u*exp(C(0, dt)); }
    f[0] = fo[0] = 0;
    char nom[96]; snprintf(nom, sizeof nom, "evol_uc%g_d%g.dat", uc, delta);
    FILE* g = fopen(nom, "w");
    fprintf(g, "# uc=%g u0=%g R=%g delta=%g\n# t  |psi(0)|/u0  N(r<4R)/N0  |psi(R)|/uc\n", uc, u0, R, delta);
    double Rn = 4*R; int iR = (int)(R/dr), in = (int)(Rn/dr);
    auto carga = [&](const vector<C>& a1, const vector<C>& a0) {   // N = ∫ Im(ψ* ψ_t) 4π r² dr = 4π ∫ Im(φ* φ_t) dr
        double s = 0; for (int i = 1; i <= in; i++) { C ph = 0.5*(a1[i] + a0[i]), pt = (a1[i] - a0[i])/dt; s += imag(conj(ph)*pt); }
        return 4*M_PI*s*dr; };
    double N0 = -carga(f, fo);
    long nsteps = (long)(tmax/dt), cada = (long)(0.5/dt);
    double minA = 1e9, maxA = 0;
    double fase_prev = 0, fase_acum = 0, fase_ini = NAN, t_ini = 0, sumN = 0, sumR = 0, sumA = 0; int cnt = 0;
    double t_prom = 0.7*tmax;                            // promedios sobre el último 30 %
    for (long n = 0; n <= nsteps; n++) {
        if (n % cada == 0) {
            double A0 = abs(f[1])/dr;                       // |ψ| cerca del centro
            double t = n*dt; if (t > 20) { minA = fmin(minA, A0/u0); maxA = fmax(maxA, A0/u0); }
            double Nt = -carga(f, fo);
            fprintf(g, "%.2f %.6f %.6f %.6f\n", t, A0/u0, Nt/N0, abs(f[iR])/(iR*dr)/uc);
            double fa = arg(f[1]); double dfa = fa - fase_prev; while (dfa > M_PI) dfa -= 2*M_PI; while (dfa < -M_PI) dfa += 2*M_PI;
            fase_acum += dfa; fase_prev = fa;
            if (t >= t_prom) {
                if (std::isnan(fase_ini)) { fase_ini = fase_acum; t_ini = t; }
                int j = 1; while (j < N-1 && abs(f[j])/(j*dr) > 0.5*A0) j++;
                sumN += Nt; sumR += j*dr; sumA += A0; cnt++;
            }
        }
        for (int i = 1; i < N-1; i++) {
            double r = i*dr; double A = abs(f[i])/r; double M = 1 - a*fmax(0.0, A - uc);
            C lap = (f[i+1] - 2.0*f[i] + f[i-1])/(dr*dr);
            C acc = lap - M*M*f[i];
            double s = sig[i]*dt;
            fn[i] = (2.0*f[i] - (1 - s)*fo[i] + dt*dt*acc)/(1 + s);
        }
        fn[0] = 0; fn[N-1] = 0;
        fo.swap(f); f.swap(fn);
    }
    fclose(g);
    printf("t > 20: |psi(0)|/u0 entre %.4f y %.4f   ->  %s\n", minA, maxA, nom);
    double tfin = nsteps*dt; double w = -(fase_acum - fase_ini)/(tfin - t_ini);
    printf("FINAL uc=%g delta=%g  N0=%.4f  N_final=%.4f  omega=%.5f  A0_medio=%.4f  R50_medio=%.3f\n",
           uc, delta, N0, sumN/cnt, w, sumA/cnt, sumR/cnt);
    return 0;
}
