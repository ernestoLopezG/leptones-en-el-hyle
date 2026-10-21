// cargafija.cpp — Envolvente con frecuencia anclada y CARGA FIJADA (prueba de la idea †, 18-10-2026)
//   u_tt = ∇²u − (M² − 1)u − σ(r)·u_t ,  M = 1 − ½(u − uc)₊      (u real: ψ = u·e^{−it}, ω anclada a 1)
//   σ(r): capa que amortigua u_t junto al borde (deja pasar la radiación sin tocar las soluciones estáticas)
//   Borde r = Rb, tres modos:
//     0  libre: capa absorbente y u = 0 muy lejos (como anclada.cpp, para comparar)
//     1  CARGA FIJA (Neumann): u_r(Rb) = −c/Rb²  ->  flujo de Gauss, Q = −4π Rb² u_r = 4π c, constante
//     2  potencial fijo (Dirichlet): u(Rb) = c/Rb
//   Uso: ./cargafija uc delta modo [tmax] [Rb]     salida: cf_uc<uc>_d<delta>_m<modo>.dat : t u(0)/u0 Q(Rb/2)/Q0 R50
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;
const double a = 0.5; double UC;
double tiro(double u0, double& R, vector<double>* X = 0, vector<double>* U = 0) {
    double uc = UC, h = 1e-4;
    auto F = [&](double x, double u, double v) { double q = 1 - a*fmax(0.0, u - uc); return -2*v/x + (q*q - 1)*u; };
    double q0 = 1 - a*(u0 - uc), F0 = (q0*q0 - 1)*u0, x = 1e-4, u = u0 + F0*x*x/6, v = F0*x/3;
    while (x < 500) {
        if (X) { X->push_back(x); U->push_back(u); }
        double k1u = v, k1v = F(x, u, v), k2u = v + .5*h*k1v, k2v = F(x + .5*h, u + .5*h*k1u, v + .5*h*k1v);
        double k3u = v + .5*h*k2v, k3v = F(x + .5*h, u + .5*h*k2u, v + .5*h*k2v), k4u = v + h*k3v, k4v = F(x + h, u + h*k3u, v + h*k3v);
        double un = u + h/6*(k1u + 2*k2u + 2*k3u + k4u), vn = v + h/6*(k1v + 2*k2v + 2*k3v + k4v);
        if (un <= uc) { double t = (u - uc)/(u - un); R = x + t*h; return R*(v + t*(vn - v)) + uc; }
        if (vn > 0) { R = NAN; return NAN; }
        u = un; v = vn; x += h;
    }
    R = NAN; return NAN;
}
int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "uso: %s uc delta modo [tmax] [Rb]\n", argv[0]); return 1; }
    UC = atof(argv[1]); double delta = atof(argv[2]); int modo = atoi(argv[3]);
    double tmax = argc > 4 ? atof(argv[4]) : 600, Rb = argc > 5 ? atof(argv[5]) : 60;
    double uc = UC, u0 = NAN, pm = NAN, pu = NAN;
    for (int i = 1; i <= 4000 && std::isnan(u0); i++) { double uu = uc + 4.0*i/4001, R; double m = tiro(uu, R);
        if (!std::isnan(m) && !std::isnan(pm) && m*pm < 0) { double lo = pu, hi = uu, ml = pm;
            for (int it = 0; it < 60; it++) { double mid = .5*(lo+hi), Rm, mm = tiro(mid, Rm); if (std::isnan(mm)) break;
                if (mm*ml < 0) hi = mid; else { lo = mid; ml = mm; } } u0 = .5*(lo+hi); }
        pm = m; pu = uu; }
    double R; vector<double> X, U; tiro(u0, R, &X, &U); double c = uc*R;
    if (modo == 0) Rb = 400;
    const double dr = 0.01, dt = 0.4*dr; int N = (int)(Rb/dr) + 1;
    auto perfil = [&](double r) { if (r >= R) return c/r; size_t j = (size_t)((r - 1e-4)/1e-4); if (j + 1 >= X.size()) j = X.size() - 2;
        double t = (r - X[j])/1e-4; return U[j] + t*(U[j+1] - U[j]); };
    vector<double> u(N), uo(N), un(N), sg(N, 0.0);
    double Ls = (modo == 0) ? 120 : 20;
    for (int i = 0; i < N; i++) { double r = i*dr; if (r > Rb - Ls) { double s = (r - (Rb - Ls))/Ls; sg[i] = 2.0*s*s; } }
    for (int i = 0; i < N; i++) u[i] = uo[i] = perfil(i ? i*dr : 1e-4)*(1 + delta);
    char nom[96]; snprintf(nom, sizeof nom, "cf_uc%g_d%g_m%d.dat", uc, delta, modo); FILE* g = fopen(nom, "w");
    fprintf(g, "# uc=%g u0=%g R=%g c=%g modo=%d Rb=%g\n# t u(0)/u0 Q(r=10)/Q0 R50\n", uc, u0, R, c, modo, Rb);
    int iq = (int)(10/dr); long ns = (long)(tmax/dt), cada = (long)(0.5/dt);
    double Amin = 1e9, Amax = -1e9;
    for (long n = 0; n <= ns; n++) {
        if (n % cada == 0) {
            double Q = -4*M_PI*(iq*dr)*(iq*dr)*(u[iq+1] - u[iq-1])/(2*dr)/(4*M_PI*c);
            int j = 1; while (j < N-1 && u[j] > 0.5*u[0]) j++;
            fprintf(g, "%.2f %.6f %.6f %.4f\n", n*dt, u[0]/u0, Q, j*dr);
            if (n*dt > 0.5*tmax) { Amin = fmin(Amin, u[0]/u0); Amax = fmax(Amax, u[0]/u0); }
        }
        for (int i = 0; i < N; i++) {
            double lap;
            if (i == 0) lap = 6*(u[1] - u[0])/(dr*dr);                         // centro: ∇²u = 3 u''(0)
            else if (i == N-1) {
                if (modo == 1) { double gh = u[N-2] - 2*dr*c/(Rb*Rb);              // u_r(Rb) = −c/Rb²
                                 lap = (gh - 2*u[i] + u[i-1])/(dr*dr) + (2/Rb)*(-c/(Rb*Rb)); }
                else { un[i] = (modo == 2) ? c/Rb : 0; continue; }
            } else { double r = i*dr; lap = (u[i+1] - 2*u[i] + u[i-1])/(dr*dr) + (u[i+1] - u[i-1])/(dr*r); }
            double M = 1 - a*fmax(0.0, u[i] - uc), s = sg[i]*dt*0.5;
            un[i] = (2*u[i] - (1 - s)*uo[i] + dt*dt*(lap - (M*M - 1)*u[i]))/(1 + s);
        }
        uo.swap(u); u.swap(un);
    }
    fclose(g);
    printf("uc=%g delta=%+g modo=%d Rb=%g : u(0)/u0 en la segunda mitad entre %.4f y %.4f\n", uc, delta, modo, Rb, Amin, Amax);
}
