// cuartico.cpp — Solitón de umbral con TÉRMINO DE CUARTO ORDEN (módulo K, abanico / Bopp–Podolsky)  († Claude, 19-10-2026)
//
// Frecuencia anclada por las compactadas (ψ = u·e^{−it}, u real). Dinámica [H]:
//     u_tt + γ·u_t = ∇²u − K·∇⁴u − U'(u) ,   U'(u) = (M² − 1)u ,  M = 1 − ½(u − uc)₊
// Energía (se conserva con γ = 0, salvo lo que se radia):
//     E = 4π ∫ [ ½u_t² + ½u_r² + (K/2)(∇²u)² + U(u) ] r² dr
// Con φ = r·u:  φ_tt = φ'' − K·φ'''' − r·U'(φ/r) − γφ_t ;  en r = 0: φ = 0 y φ'' = 0 (extensión impar).
// Escalas (Derrick con K): gradiente ∝ R, término K ∝ 1/R, potencial ∝ R³. El término K cuesta al encoger:
//   puede crear un mínimo local de energía (estado estable frente a empujones pequeños).
//
// Uso:  ./cuartico uc K gamma delta tmax [entrada.dat]
//   sin entrada: arranca de la solución estática de umbral (umbral.cpp) multiplicada por (1+delta)
//   con entrada: arranca del perfil guardado (r u) multiplicado por (1+delta)
// Salidas: cu_<etiqueta>.dat : t  u(0)  R50  E   ;  final_<etiqueta>.dat : perfil final (r u)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
using namespace std;
double uc;
double Uprime(double u) { double M = 1 - 0.5*fmax(0.0, u - uc); return (M*M - 1)*u; }
double U(double u) { if (u <= uc) return 0; double Y = u - uc; return -Y*Y*Y/3 - uc*Y*Y/2 + Y*Y*Y*Y/16 + uc*Y*Y*Y/12; }
double tiro(double u0, double& R, vector<double>* X = 0, vector<double>* Uv = 0) {
    double h = 1e-4; auto F = [&](double x, double u, double v) { return -2*v/x + Uprime(u); };
    double F0 = Uprime(u0), x = 1e-4, u = u0 + F0*x*x/6, v = F0*x/3;
    while (x < 500) { if (X) { X->push_back(x); Uv->push_back(u); }
        double k1u = v, k1v = F(x,u,v), k2u = v+.5*h*k1v, k2v = F(x+.5*h,u+.5*h*k1u,v+.5*h*k1v);
        double k3u = v+.5*h*k2v, k3v = F(x+.5*h,u+.5*h*k2u,v+.5*h*k2v), k4u = v+h*k3v, k4v = F(x+h,u+h*k3u,v+h*k3v);
        double un = u + h/6*(k1u+2*k2u+2*k3u+k4u), vn = v + h/6*(k1v+2*k2v+2*k3v+k4v);
        if (un <= uc) { double t = (u-uc)/(u-un); R = x + t*h; return R*(v + t*(vn-v)) + uc; }
        if (vn > 0) { R = NAN; return NAN; } u = un; v = vn; x += h; }
    R = NAN; return NAN; }
int main(int argc, char** argv) {
    uc = atof(argv[1]); double K = atof(argv[2]), gam = atof(argv[3]), delta = atof(argv[4]), tmax = atof(argv[5]);
    const double h = 0.05, L = 150, Ls = 40; int N = (int)(L/h) + 1;
    double kmax = M_PI/h, wmax = sqrt(kmax*kmax + K*pow(kmax,4) + 10), dt = 0.5/wmax;
    vector<double> r(N), ui(N);
    for (int i = 0; i < N; i++) r[i] = i*h;
    if (argc > 6) {                                    // perfil de entrada
        FILE* f = fopen(argv[6], "r"); vector<double> rr, uu; double a, b; char line[256];
        while (fgets(line, sizeof line, f)) if (line[0] != '#' && sscanf(line, "%lf %lf", &a, &b) == 2) { rr.push_back(a); uu.push_back(b); }
        fclose(f);
        for (int i = 0; i < N; i++) { size_t j = 0; while (j + 1 < rr.size() && rr[j+1] < r[i]) j++;
            ui[i] = (j + 1 < rr.size()) ? uu[j] + (uu[j+1]-uu[j])*(r[i]-rr[j])/(rr[j+1]-rr[j]) : uu.back()*rr.back()/fmax(r[i],1e-9); }
    } else {
        double u0 = NAN, pm = NAN, pu = NAN;
        for (int i = 1; i <= 4000 && std::isnan(u0); i++) { double uu = uc + 4.0*i/4001, R; double m = tiro(uu, R);
            if (!std::isnan(m) && !std::isnan(pm) && m*pm < 0) { double lo = pu, hi = uu, ml = pm;
                for (int it = 0; it < 60; it++) { double mid = .5*(lo+hi), Rm, mm = tiro(mid, Rm); if (std::isnan(mm)) break;
                    if (mm*ml < 0) hi = mid; else { lo = mid; ml = mm; } } u0 = .5*(lo+hi); }
            pm = m; pu = uu; }
        double R; vector<double> X, Uv; tiro(u0, R, &X, &Uv); double c = uc*R;
        for (int i = 0; i < N; i++) { double x = fmax(r[i], 1e-4);
            if (x >= R) ui[i] = c/x; else { size_t j = (size_t)((x - 1e-4)/1e-4); if (j >= Uv.size()) j = Uv.size()-1; ui[i] = Uv[j]; } }
    }
    vector<double> p(N), po(N), pn(N), sg(N, 0);
    for (int i = 0; i < N; i++) { p[i] = po[i] = r[i]*ui[i]*(1 + delta); if (r[i] > L - Ls) sg[i] = 0.5*pow((r[i]-(L-Ls))/Ls, 2); }
    char et[128]; snprintf(et, sizeof et, "uc%g_K%g_g%g_d%g%s", uc, K, gam, delta, argc > 6 ? "_rel" : "");
    string nom = string("cu_") + et + ".dat"; FILE* g = fopen(nom.c_str(), "w");
    fprintf(g, "# uc=%g K=%g gamma=%g delta=%g dt=%g\n# t u(0) R50 E\n", uc, K, gam, delta, dt);
    auto P = [&](const vector<double>& a, int i) { return i >= 0 ? (i < N ? a[i] : 0.0) : -a[-i]; };   // impar en r=0
    long ns = (long)(tmax/dt), cada = (long)(1.0/dt);
    for (long n = 0; n <= ns; n++) {
        if (n % cada == 0) {
            double u0 = 2*p[1]/h - p[2]/(2*h);                         // u(0) ≈ φ'(0)
            int j = 1; while (j < N-1 && p[j]/r[j] > 0.5*u0) j++;
            double E = 0;
            for (int i = 1; i < N-1; i++) { double ri = r[i], u = p[i]/ri, ut = (p[i]-po[i])/(dt*ri);
                double ur = ((p[i+1]/r[i+1]) - (i > 1 ? p[i-1]/r[i-1] : u0))/(2*h);
                double lap = (P(p,i+1) - 2*p[i] + P(p,i-1))/(h*h)/ri;
                if (sg[i] > 0) continue;
                E += 4*M_PI*(0.5*ut*ut + 0.5*ur*ur + 0.5*K*lap*lap + U(u))*ri*ri*h; }
            fprintf(g, "%.2f %.6f %.3f %.5f\n", n*dt, u0, j*h, E); fflush(g);
        }
        for (int i = 1; i < N-1; i++) {
            double d2 = (P(p,i+1) - 2*p[i] + P(p,i-1))/(h*h);
            double d4 = (P(p,i+2) - 4*P(p,i+1) + 6*p[i] - 4*P(p,i-1) + P(p,i-2))/(h*h*h*h);
            double acc = d2 - K*d4 - r[i]*Uprime(p[i]/r[i]);
            double s = (gam + sg[i])*dt*0.5;
            pn[i] = (2*p[i] - (1 - s)*po[i] + dt*dt*acc)/(1 + s);
        }
        pn[0] = 0; pn[N-1] = 0; po.swap(p); p.swap(pn);
    }
    fclose(g);
    string fin = string("final_") + et + ".dat"; FILE* f = fopen(fin.c_str(), "w"); fprintf(f, "# r u\n");
    fprintf(f, "%g %g\n", 0.0, 2*p[1]/h - p[2]/(2*h));
    for (int i = 1; i < N; i++) fprintf(f, "%g %g\n", r[i], p[i]/r[i]); fclose(f);
    printf("%s : dt=%.2e  ", et, dt);
    FILE* q = fopen(nom.c_str(), "r"); char line[256]; double t, a, R, E, tl = 0, al = 0, Rl = 0, El = 0;
    while (fgets(line, sizeof line, q)) if (line[0] != '#' && sscanf(line, "%lf %lf %lf %lf", &t, &a, &R, &E) == 4) { tl = t; al = a; Rl = R; El = E; }
    fclose(q); printf("final t=%.0f u(0)=%.4f R50=%.3f E=%.4f\n", tl, al, Rl, El);
}
