// hopfion_compacto.cpp — ¿llena el nudo una dimensión compacta, y ajusta entonces su tamaño en las extendidas?
// (Claude, 4-11-2026)
//
// Pregunta (Ernesto): «¿qué limita el tamaño del nudo? En principio nada, pero la forma sí:
//   me expando hasta tener el mismo tamaño en todas las dimensiones».
// Precedente: un hopfión en una 3-esfera de radio R la llena entera si R < √2 y se localiza si R > √2
//   (R. S. Ward, «Hopf solitons on S³ and R³», arXiv:hep-th/9811176).
//
// Modelo de juguete: el mismo hopfión de Faddeev–Skyrme que hopfion.cpp,
//     E = ∫ [ ∂_i n·∂_i n + ½ F_ij F_ij ] d³x ,  F_ij = n·(∂_i n × ∂_j n),
// pero en R² × S¹: x, y extendidas (caja de lado L, borde fijo n = (0,0,1)) y z COMPACTA, periódica de periodo Lc.
// z hace el papel de una dirección compacta (η). No es el hyle completo: solo prueba el mecanismo.
//
// Método: continuación. Se arranca con el hopfión de carga 1 (tamaño a) en un periodo grande Lc₁,
// se relaja, y luego se va encogiendo Lc paso a paso (se reescala z: mismos valores en la malla,
// menor espaciado hz = Lc/Nz) y se vuelve a relajar. La carga de Hopf se conserva en la deformación.
//
// Uso:   ./hopfion_compacto N Nz L a itmax Lc1 Lc2 ... [-t etiqueta]
//   N    : puntos en x e y        Nz : puntos en z (periodo)       L : lado de la caja en x, y
//   a    : tamaño del hopfión inicial          itmax : iteraciones máximas por etapa
//   Lc1 > Lc2 > ... : periodos de z que se recorren (de mayor a menor)
//   Variables de entorno: DTAU (paso, por defecto 0.02·h²), TOL (parada, por defecto 1e-6 de cambio relativo de E
//   por cada 500 iteraciones), ITMIN (mínimo por etapa, por defecto 2000).
// Salidas:
//   compacto_<etiqueta>.dat : Lc  iter  E2  E4  E/(32√2π²)  R_xy  sigma_z  Rbar  min/max(z)
//   n_<etiqueta>_Lc<Lc>.bin : estado final de cada etapa (para hopfq_compacto.py y figuras_compacto.py)
// Medidas (pesando con w = 1 − n_z, como medir_tamano.py):
//   R_xy    : radio cuadrático medio en el plano extendido, √(Σ w (x²+y²) / Σ w)
//   Rbar    : concentración en z (media circular): 1 = localizado en z, 0 = repartido uniformemente
//   sigma_z : anchura en z (desviación circular, (Lc/2π)·√(−2 ln Rbar)); comparable con R_xy/√2 si el nudo es redondo
//   min/max : perfil Σ_xy w a lo largo de z, mínimo/máximo; ≈ 1 si el nudo llena toda la dirección compacta
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>
using namespace std;
struct V { double x, y, z; };
inline V operator+(V a, V b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline V operator-(V a, V b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline V operator*(double s, V a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(V a, V b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline V cross(V a, V b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }

int N, Nz;
inline size_t I(int i, int j, int k) { k = (k + Nz) % Nz; return ((size_t)i*N + j)*Nz + k; }

int main(int argc, char** argv) {
    if (argc < 7) { fprintf(stderr, "uso: ./hopfion_compacto N Nz L a itmax Lc1 Lc2 ... [-t etiqueta]\n"); return 1; }
    N = atoi(argv[1]); Nz = atoi(argv[2]); double L = atof(argv[3]), a = atof(argv[4]); int itmax = atoi(argv[5]);
    vector<double> Lcs; string tag = "c";
    for (int i = 6; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i+1 < argc) { tag = argv[++i]; continue; }
        Lcs.push_back(atof(argv[i]));
    }
    const double hx = L/(N-1);
    const double dtf = getenv("DTAU") ? atof(getenv("DTAU")) : 0.02;
    const double tol = getenv("TOL") ? atof(getenv("TOL")) : 1e-6;
    const int itmin = getenv("ITMIN") ? atoi(getenv("ITMIN")) : 2000;
    const double norma = 32*sqrt(2.0)*M_PI*M_PI;
    size_t NT = (size_t)N*N*Nz;
    vector<V> n(NT), G(NT), d[3], W[3];
    for (int p = 0; p < 3; p++) { d[p].resize(NT); W[p].resize(NT); }

    // Estado inicial: hopfión de Hopf de carga 1, tamaño a, centrado, en el primer periodo Lc1.
    double hz = Lcs[0]/Nz;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < Nz; k++) {
        double x = (-L/2 + i*hx)/a, y = (-L/2 + j*hx)/a, z = (-Lcs[0]/2 + (k+0.5)*hz)/a, r2 = x*x + y*y + z*z;
        double u1 = 2*x/(1+r2), u2 = 2*y/(1+r2), u3 = 2*z/(1+r2), u4 = (r2-1)/(1+r2);
        double z1r = u1, z1i = u2, z0r = u4, z0i = u3;
        double re = z1r*z0r + z1i*z0i, im = z1i*z0r - z1r*z0i, m1 = z1r*z1r + z1i*z1i, m0 = z0r*z0r + z0i*z0i;
        n[I(i,j,k)] = {2*re/(m0+m1), 2*im/(m0+m1), (m0-m1)/(m0+m1)};
        if (i == 0 || j == 0 || i == N-1 || j == N-1) n[I(i,j,k)] = {0, 0, 1};
    }

    string nomdat = "compacto_" + tag + ".dat";
    FILE* f = fopen(nomdat.c_str(), "w");
    fprintf(f, "# N=%d Nz=%d L=%g a=%g hx=%g DTAU=%g TOL=%g\n# Lc iter E2 E4 E/(32√2π²) R_xy sigma_z Rbar min/max(z)\n",
            N, Nz, L, a, hx, dtf, tol);
    printf("hopfión en R² × S¹: N=%d, Nz=%d, L=%g, a=%g, %zu etapas\n", N, Nz, L, a, Lcs.size());

    for (size_t etapa = 0; etapa < Lcs.size(); etapa++) {
        double Lc = Lcs[etapa]; hz = Lc/Nz;
        double hmin = min(hx, hz), dtau = dtf*hmin*hmin;
        double E2 = 0, E4 = 0, Eant = 1e300; int it;
        for (it = 0; it <= itmax; it++) {
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 0; k < Nz; k++) {
                size_t c = I(i,j,k);
                d[0][c] = (0.5/hx)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]);
                d[1][c] = (0.5/hx)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]);
                d[2][c] = (0.5/hz)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
            }
            E2 = 0; E4 = 0;
            #pragma omp parallel for collapse(2) reduction(+:E2,E4)
            for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 0; k < Nz; k++) {
                size_t c = I(i,j,k); V nn = n[c]; double F[3][3];
                for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) F[p][q] = dot(nn, cross(d[p][c], d[q][c]));
                for (int p = 0; p < 3; p++) { V s = {0,0,0};
                    for (int q = 0; q < 3; q++) s = s + F[p][q]*cross(d[q][c], nn);
                    W[p][c] = s; E2 += dot(d[p][c], d[p][c]); }
                E4 += F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2];
            }
            E2 *= hx*hx*hz; E4 *= hx*hx*hz;
            if (it % 500 == 0) {
                double E = E2 + E4;
                if (it >= itmin && fabs(Eant - E)/E < tol) break;
                Eant = E;
                if (it % 5000 == 0) printf("  Lc %6.3f  it %6d  E/(32√2π²) %.5f  E2/E4 %.3f\n", Lc, it, E/norma, E2/E4);
                fflush(stdout);
            }
            if (it == itmax) break;
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 0; k < Nz; k++) {
                size_t c = I(i,j,k);
                V lap = (1/(hx*hx))*(n[I(i+1,j,k)] + n[I(i-1,j,k)] + n[I(i,j+1,k)] + n[I(i,j-1,k)] - 4.0*n[c])
                      + (1/(hz*hz))*(n[I(i,j,k+1)] + n[I(i,j,k-1)] - 2.0*n[c]);
                V divW = (0.5/hx)*(W[0][I(i+1,j,k)] - W[0][I(i-1,j,k)] + W[1][I(i,j+1,k)] - W[1][I(i,j-1,k)])
                       + (0.5/hz)*(W[2][I(i,j,k+1)] - W[2][I(i,j,k-1)]);
                V g = -2.0*lap - 2.0*divW; g = g - dot(g, n[c])*n[c];
                G[c] = g;
            }
            #pragma omp parallel for collapse(2)
            for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 0; k < Nz; k++) {
                size_t c = I(i,j,k); V m = n[c] - dtau*G[c]; n[c] = (1/sqrt(dot(m,m)))*m; }
        }
        // medidas
        double sw = 0, sr2 = 0; vector<double> wz(Nz, 0.0);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < Nz; k++) {
            double x = -L/2 + i*hx, y = -L/2 + j*hx, w = 1 - n[I(i,j,k)].z;
            sw += w; sr2 += w*(x*x + y*y); wz[k] += w;
        }
        double C = 0, S = 0;
        for (int k = 0; k < Nz; k++) { double th = 2*M_PI*(k+0.5)/Nz; C += wz[k]*cos(th); S += wz[k]*sin(th); }
        double Rbar = sqrt(C*C + S*S)/sw;
        double sz = Rbar > 1e-12 ? Lc/(2*M_PI)*sqrt(-2*log(Rbar)) : INFINITY;
        double mm = *min_element(wz.begin(), wz.end()) / *max_element(wz.begin(), wz.end());
        double Rxy = sqrt(sr2/sw);
        fprintf(f, "%.4f %d %.5f %.5f %.5f %.4f %.4f %.4f %.4f\n", Lc, it, E2, E4, (E2+E4)/norma, Rxy, sz, Rbar, mm); fflush(f);
        printf("Lc %6.3f  iter %6d  E/(32√2π²) %.5f  R_xy %.3f  sigma_z %.3f  Rbar %.3f  min/max %.3f\n",
               Lc, it, (E2+E4)/norma, Rxy, sz, Rbar, mm);
        char nb[128]; snprintf(nb, sizeof nb, "n_%s_Lc%.3f.bin", tag.c_str(), Lc);
        FILE* b = fopen(nb, "wb");
        fwrite(&N, sizeof N, 1, b); fwrite(&Nz, sizeof Nz, 1, b); fwrite(&L, sizeof L, 1, b); fwrite(&Lc, sizeof Lc, 1, b);
        fwrite(n.data(), sizeof(V), NT, b); fclose(b);
    }
    fclose(f);
    printf("-> %s\n", nomdat.c_str());
}
