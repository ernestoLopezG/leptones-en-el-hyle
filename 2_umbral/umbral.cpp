// umbral.cpp — Leptón en reposo con umbral inferior tipo Bingham (opción (a) de Ernesto, 18-10-2026)  († Claude)
//
//   El hyle no responde mientras u < uc (umbral de fluencia). Por encima responde como antes:
//       u'' + (2/x) u' = [ q² − 1 ] u ,   q = 1 − a·(u − uc)₊ ,   a = 1/2 ,  eps = 0, sin autocampo (gw = 0)
//   Fuera (u < uc): ∇²u = 0  ->  u = c/x  (Coulomb puro, lineal: orbitales y paquetes viven aquí).
//   Empalme en R (u = uc) con u y u' continuas:  c = uc·R  y  u'(R) = −uc/R   <- condición que FIJA u0.
//   Techo (meseta de Klein, q = −1):  u = uc + 4.
//   Carga ρ = −∇²u = (1 − q²)u, solo dentro (u > uc); Q = 4π c = 4π uc R  (Gauss).
// Uso:  ./umbral uc1 uc2 ...      -> tabla; perfil_uc<uc>.dat (x u rho) de la primera solución
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;
const double a = 0.5, dx = 1e-4;

// integra desde el centro hasta u = uc; devuelve mismatch m = R·u'(R) + uc (0 = empalme) ; R por referencia
double tiro(double u0, double uc, double& R, vector<double>* X = 0, vector<double>* U = 0, vector<double>* RH = 0) {
    auto F = [&](double x, double u, double v) { double q = 1 - a*fmax(0.0, u - uc); return -2*v/x + (q*q - 1)*u; };
    double q0 = 1 - a*(u0 - uc), F0 = (q0*q0 - 1)*u0;
    double x = 1e-3, u = u0 + F0*x*x/6, v = F0*x/3;
    while (x < 200) {
        if (X) { double q = 1 - a*fmax(0.0, u - uc); X->push_back(x); U->push_back(u); RH->push_back((1 - q*q)*u); }
        double k1u = v, k1v = F(x, u, v);
        double k2u = v + .5*dx*k1v, k2v = F(x + .5*dx, u + .5*dx*k1u, v + .5*dx*k1v);
        double k3u = v + .5*dx*k2v, k3v = F(x + .5*dx, u + .5*dx*k2u, v + .5*dx*k2v);
        double k4u = v + dx*k3v,    k4v = F(x + dx, u + dx*k3u, v + dx*k3v);
        double un = u + dx/6*(k1u + 2*k2u + 2*k3u + k4u), vn = v + dx/6*(k1v + 2*k2v + 2*k3v + k4v);
        if (un <= uc) { double t = (u - uc)/(u - un); R = x + t*dx; double vR = v + t*(vn - v); return R*vR + uc; }
        if (vn > 0) { R = NAN; return NAN; }      // vuelve a subir: no llega al umbral
        u = un; v = vn; x += dx;
    }
    R = NAN; return NAN;
}

int main(int argc, char** argv) {
    printf("#  uc      u0        R       c=uc·R   Q=4πc     R50      rq(carga)   [x en ħ/mc]\n");
    for (int k = 1; k < argc; k++) {
        double uc = atof(argv[k]); int n = 0; bool primero = true;
        double prevm = NAN, prevu = NAN;
        for (int i = 1; i <= 4000; i++) {                       // barrido de u0 en (uc, uc+4)
            double u0 = uc + 4.0*i/4001, R; double m = tiro(u0, uc, R);
            if (!std::isnan(m) && !std::isnan(prevm) && m*prevm < 0) {
                double lo = prevu, hi = u0, mlo = prevm;
                for (int it = 0; it < 60; it++) { double mid = .5*(lo+hi), Rm; double mm = tiro(mid, uc, Rm);
                    if (std::isnan(mm)) break; if (mm*mlo < 0) hi = mid; else { lo = mid; mlo = mm; } }
                double u0s = .5*(lo+hi), Rs; tiro(u0s, uc, Rs);
                vector<double> X, U, RH; tiro(u0s, uc, Rs, &X, &U, &RH);
                // perfil completo con cola c/x, R50 y radio de carga
                double c = uc*Rs, R50 = NAN, Q = 0, r2 = 0;
                for (size_t j = 0; j < X.size(); j++) { if (std::isnan(R50) && U[j] <= .5*u0s) R50 = X[j];
                    Q += RH[j]*4*M_PI*X[j]*X[j]*dx; r2 += RH[j]*X[j]*X[j]*4*M_PI*X[j]*X[j]*dx; }
                if (std::isnan(R50)) R50 = c/(.5*u0s);           // la mitad cae en la cola
                printf("%6.3f %9.5f %8.4f %9.4f %9.4f %8.4f %9.4f   (Q dentro = %.4f)\n", uc, u0s, Rs, c, 4*M_PI*c, R50, sqrt(r2/Q), Q);
                if (primero) { char nom[64]; snprintf(nom, sizeof nom, "perfil_uc%g.dat", uc); FILE* f = fopen(nom, "w");
                    for (size_t j = 0; j < X.size(); j += 50) fprintf(f, "%g %g %g\n", X[j], U[j], RH[j]);
                    for (double x = Rs; x < 6*Rs + 20; x += 0.05) fprintf(f, "%g %g 0\n", x, c/x);
                    fclose(f); primero = false; }
                n++;
            }
            prevm = m; prevu = u0;
        }
        if (!n) printf("%6.3f   sin empalme en (uc, uc+4)\n", uc);
    }
}
