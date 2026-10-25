// giro_anc.cpp — como giro.cpp, con opción ANCLADA=1 (frecuencia anclada por las compactadas). Envolvente que GIRA alrededor del eje z (simetría de eje, coordenadas ρ, z)  († Claude, 18-10-2026)
//
//   ψ(ρ,z,φ,t) = w(ρ,z,t)·e^{i m φ}     (m = número de vueltas del giro; m = 0 sin giro)
//   Ecuación [H] (la de evolucion.cpp, frecuencia libre, con umbral inferior uc):
//       w_tt = w_ρρ + w_ρ/ρ − m² w/ρ² + w_zz − M² w ,   M = 1 − ½(|w| − uc)₊
//   Con giro se conservan la carga de fase N = ∫ Im(ψ* ψ_t) d³x... y el momento angular J = m·N.
//   Preguntas: ¿forma anillo?, ¿aguanta?, ¿a qué frecuencia gira la fase (ω = 1 sin autoenlace; ω < 1 autoenlace)?
//
//   Arranque: w = A·(ρ/w0)^m·exp(−(ρ² + z²)/w0²)·e^{−it}   (un bulto que ya gira)
//   Uso: [ANCLADA=1] ./giro_anc m A w0 tmax uc          salida: giro_m<m>_A<A>.dat :
//        t  |w|max  ρ_max  z-ancho(media altura)  N/N0  ω (fase en el máximo)
//        y perfiles |w|(ρ,z) en t = 0, tmax/4, tmax/2, tmax  (giro_m<m>_A<A>_t<t>.dat)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <complex>
#include <omp.h>
using namespace std;
bool ANC = getenv("ANCLADA") != nullptr;   // ANCLADA=1: ψ = U·e^{imφ}·e^{−it}, U real, frecuencia anclada (parte real)
typedef complex<double> C;
int main(int argc, char** argv) {
    int m = atoi(argv[1]); double A = atof(argv[2]), w0 = atof(argv[3]), tmax = atof(argv[4]), uc = atof(argv[5]);
    const double h = 0.1, dt = 0.04, Lr = 50, Lz = 50, Ls = 15;
    int NR = (int)(Lr/h) + 1, NZ = (int)(2*Lz/h) + 1;
    auto id = [&](int i, int k) { return (size_t)i*NZ + k; };
    vector<C> w(NR*NZ), wo(NR*NZ), wn(NR*NZ); vector<double> sg(NR*NZ, 0);
    for (int i = 0; i < NR; i++) for (int k = 0; k < NZ; k++) {
        double r = i*h, z = -Lz + k*h, d = fmax(r - (Lr - Ls), fabs(z) - (Lz - Ls));
        if (d > 0) sg[id(i,k)] = 2.0*pow(d/Ls, 2);
        double a = A*pow(r/w0, m)*exp(-(r*r + z*z)/(w0*w0));
        w[id(i,k)] = a; wo[id(i,k)] = ANC ? C(a) : a*exp(C(0, dt));         // e^{−it}: paso atrás
    }
    auto carga = [&]() { double s = 0;
        #pragma omp parallel for reduction(+:s)
        for (int i = 1; i < NR-1; i++) for (int k = 1; k < NZ-1; k++) { size_t j = id(i,k);
            if (sg[j] > 0) continue; C wm = 0.5*(w[j] + wo[j]), wt = (w[j] - wo[j])/dt;
            s += imag(conj(wm)*wt)*2*M_PI*i*h*h*h; }
        return -s; };
    char nom[96]; snprintf(nom, sizeof nom, ANC ? "anc_m%d_A%g.dat" : "giro_m%d_A%g.dat", m, A); FILE* g = fopen(nom, "w");
    fprintf(g, "# m=%d A=%g w0=%g uc=%g\n# t |w|max rho_max z_ancho N/N0 omega\n", m, A, w0, uc);
    double N0 = carga(), fase_prev = 0; long ns = (long)(tmax/dt), cada = (long)(1.0/dt);
    for (long n = 0; n <= ns; n++) {
        double t = n*dt;
        if (n % cada == 0) {
            int im = 0, km = 0; double amax = 0;
            for (int i = 0; i < NR; i++) for (int k = 0; k < NZ; k++) { double a = abs(w[id(i,k)]); if (a > amax) { amax = a; im = i; km = k; } }
            int k1 = km, k2 = km; while (k1 > 0 && abs(w[id(im,k1)]) > .5*amax) k1--; while (k2 < NZ-1 && abs(w[id(im,k2)]) > .5*amax) k2++;
            double fa = arg(w[id(im,km)]), om = -(fa - fase_prev); while (om > M_PI) om -= 2*M_PI; while (om < -M_PI) om += 2*M_PI;
            fase_prev = fa;
            fprintf(g, "%.1f %.5f %.2f %.2f %.5f %.4f\n", t, amax, im*h, (k2 - k1)*h, carga()/N0, om/(cada*dt)); fflush(g);
        }
        for (double tt : {0.0, tmax/4, tmax/2, tmax}) if (fabs(t - tt) < 0.5*dt) {
            char pn[96]; snprintf(pn, sizeof pn, ANC ? "anc_m%d_A%g_t%g.dat" : "giro_m%d_A%g_t%g.dat", m, A, tt); FILE* p = fopen(pn, "w");
            for (int i = 0; i < NR; i += 2) { for (int k = 0; k < NZ; k += 2) fprintf(p, "%.4f ", abs(w[id(i,k)])); fprintf(p, "\n"); }
            fclose(p); }
        #pragma omp parallel for
        for (int i = 0; i < NR; i++) for (int k = 0; k < NZ; k++) {
            size_t j = id(i,k);
            if (i == NR-1 || k == 0 || k == NZ-1) { wn[j] = 0; continue; }
            C lap;
            if (i == 0) { if (m != 0) { wn[j] = 0; continue; }
                lap = 4.0*(w[id(1,k)] - w[j])/(h*h) + (w[id(0,k+1)] - 2.0*w[j] + w[id(0,k-1)])/(h*h); }
            else { double r = i*h;
                lap = (w[id(i+1,k)] - 2.0*w[j] + w[id(i-1,k)])/(h*h) + (w[id(i+1,k)] - w[id(i-1,k)])/(2*h*r)
                    - double(m*m)*w[j]/(r*r) + (w[id(i,k+1)] - 2.0*w[j] + w[id(i,k-1)])/(h*h); }
            double M = 1 - 0.5*fmax(0.0, abs(w[j]) - uc), s = sg[j]*dt*0.5;
            wn[j] = (2.0*w[j] - (1 - s)*wo[j] + dt*dt*(lap - (ANC ? M*M - 1.0 : M*M)*w[j]))/(1 + s);
        }
        wo.swap(w); w.swap(wn);
    }
    fclose(g);
}
