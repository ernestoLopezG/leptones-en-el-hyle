// derrick.cpp — ¿Es inestable la solución de reposo2 SIN umbral inferior? (frecuencia anclada, u real)  († Claude, 18-10)
//   Estática (reposo2, eps = 0, gw = 0): u'' + (2/x)u' = [(1 − u/2)² − 1]u, empalmada con c/x en x_b (u + x·u' = 0)
//   Dinámica: u_tt = ∇²u − [(1 − u/2)² − 1]u   (la misma de anclada.cpp con uc = 0)
//   Uso: ./derrick u0 delta tmax    -> imprime u(0)/u0 cada 5 unidades de tiempo
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;
int main(int argc, char** argv) {
    double u0 = atof(argv[1]), delta = atof(argv[2]), tmax = atof(argv[3]);
    double h = 1e-4, x = 1e-4, F0 = (pow(1 - u0/2, 2) - 1)*u0, u = u0 + F0*x*x/6, v = F0*x/3, xb = NAN, c = 0;
    vector<double> X, U;
    while (x < 100) { X.push_back(x); U.push_back(u);
        auto F = [&](double xx, double uu, double vv) { return -2*vv/xx + (pow(1 - uu/2, 2) - 1)*uu; };
        double k1u = v, k1v = F(x, u, v), k2u = v + .5*h*k1v, k2v = F(x + .5*h, u + .5*h*k1u, v + .5*h*k1v);
        double k3u = v + .5*h*k2v, k3v = F(x + .5*h, u + .5*h*k2u, v + .5*h*k2v), k4u = v + h*k3v, k4v = F(x + h, u + h*k3u, v + h*k3v);
        u += h/6*(k1u + 2*k2u + 2*k3u + k4u); v += h/6*(k1v + 2*k2v + 2*k3v + k4v); x += h;
        if (u + x*v <= 0) { xb = x; c = x*u; break; } }
    printf("# u0=%g  x_b=%.4f  c=%.4f  delta=%g\n", u0, xb, c, delta);
    const double dr = 0.01, dt = 0.4*dr, L = 400, Ls = 120; int N = (int)(L/dr) + 1;
    vector<double> w(N), wo(N), wn(N), sg(N, 0);
    for (int i = 0; i < N; i++) { double r = i ? i*dr : 1e-4, val;
        if (r >= xb) val = c/r; else { size_t j = (size_t)((r - 1e-4)/h); if (j >= U.size()) j = U.size() - 1; val = U[j]; }
        w[i] = wo[i] = val*(1 + delta);
        if (r > L - Ls) { double s = (r - (L - Ls))/Ls; sg[i] = 2*s*s; } }
    long ns = (long)(tmax/dt), cada = (long)(5/dt);
    for (long n = 0; n <= ns; n++) {
        if (n % cada == 0) printf("%5.0f %8.4f\n", n*dt, w[0]/u0);
        for (int i = 0; i < N - 1; i++) {
            double lap = i ? (w[i+1] - 2*w[i] + w[i-1])/(dr*dr) + (w[i+1] - w[i-1])/(dr*i*dr) : 6*(w[1] - w[0])/(dr*dr);
            double s = sg[i]*dt*.5, f = (pow(1 - w[i]/2, 2) - 1)*w[i];
            wn[i] = (2*w[i] - (1 - s)*wo[i] + dt*dt*(lap - f))/(1 + s); }
        wn[N-1] = 0; wo.swap(w); w.swap(wn); }
}
