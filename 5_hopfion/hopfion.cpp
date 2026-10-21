// hopfion.cpp — Paso 1: hopfión de referencia (modelo de Faddeev–Skyrme) para validar el método  († Claude, 19-10-2026)
//
// Campo: director unitario n(x,y,z) ∈ S² (en el hyle: la orientación del director).
// Energía (Faddeev–Skyrme; convención de Battye y Sutcliffe):
//     E = ∫ [ ∂_i n·∂_i n  +  ½ F_ij F_ij ] d³x ,   F_ij = n·(∂_i n × ∂_j n)
//   2.º término = término de cuarto orden (tipo Skyrme / nuestro K): frena el colapso.
//   Referencia: en la normalización E/(32√2π²), E(Q=1) ≈ 1,21 (Sutcliffe, Proc. R. Soc. A 463, 3001 (2007), arXiv:0705.1468).
//   (El valor 504 de Battye y Sutcliffe, PRL 81, 4798 (1998), malla 100³, resultó ser una subestimación.)
//   Comprobación independiente de la norma: en el mínimo, E2 = E4 (teorema del virial / Derrick).
// Método: flujo de gradiente proyectado sobre la esfera; n fijo = (0,0,1) en el borde.
// Arranque: mapa de Hopf de carga 1, tamaño a.
// Uso: ./hopfion N L a iter [continuar.bin]   (malla N³ en [−L/2, L/2]³; con .bin continúa desde ese estado,
//      que debe tener la misma N y L). Cada 2000 iteraciones guarda n_final_N<N>.bin (se puede cortar y seguir).
//    salida: hopf_N<N>.dat (iter E2 E4 E/(32√2π²) E2/E4)
//      y n_final_N<N>.bin para calcular Q con hopfq.py
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <omp.h>
using namespace std;
struct V { double x, y, z; };
inline V operator+(V a, V b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline V operator-(V a, V b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline V operator*(double s, V a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(V a, V b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline V cross(V a, V b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
int N; double h;
inline size_t I(int i, int j, int k) { return ((size_t)i*N + j)*N + k; }
int main(int argc, char** argv) {
    N = atoi(argv[1]); double L = atof(argv[2]), a = atof(argv[3]); int iters = atoi(argv[4]);
    h = L/(N-1);
    vector<V> n(N*N*N), G(N*N*N); vector<V> d[3]; for (auto& v : d) v.resize(N*N*N);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        double x = (-L/2 + i*h)/a, y = (-L/2 + j*h)/a, z = (-L/2 + k*h)/a, r2 = x*x + y*y + z*z;
        double u1 = 2*x/(1+r2), u2 = 2*y/(1+r2), u3 = 2*z/(1+r2), u4 = (r2-1)/(1+r2);   // R³ -> S³
        double z1r = u1, z1i = u2, z0r = u4, z0i = u3;                                 // Z1 = u1+iu2, Z0 = u4+iu3
        double re = z1r*z0r + z1i*z0i, im = z1i*z0r - z1r*z0i, m1 = z1r*z1r + z1i*z1i, m0 = z0r*z0r + z0i*z0i;
        n[I(i,j,k)] = {2*re/(m0+m1), 2*im/(m0+m1), (m0-m1)/(m0+m1)};
        if (i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1) n[I(i,j,k)] = {0, 0, 1};
    }
    if (argc > 5) { FILE* r = fopen(argv[5], "rb"); int N2; double L2;
        if (!r || fread(&N2, sizeof N2, 1, r) != 1 || fread(&L2, sizeof L2, 1, r) != 1 || N2 != N || fabs(L2 - L) > 1e-9) {
            fprintf(stderr, "el fichero de continuación no coincide en N o L\n"); return 1; }
        if (fread(n.data(), sizeof(V), n.size(), r) != n.size()) { fprintf(stderr, "fichero incompleto\n"); return 1; }
        fclose(r); printf("continuando desde %s\n", argv[5]); }
    auto guarda = [&]() { char nb[64]; snprintf(nb, sizeof nb, "n_final_N%d.bin", N); FILE* b = fopen(nb, "wb");
        fwrite(&N, sizeof N, 1, b); fwrite(&L, sizeof L, 1, b); fwrite(n.data(), sizeof(V), n.size(), b); fclose(b); };
    char nom[64]; snprintf(nom, sizeof nom, "hopf_N%d.dat", N); FILE* f = fopen(nom, argc > 5 ? "a" : "w");
    fprintf(f, "# N=%d L=%g a=%g h=%g\n# iter E2 E4 E/(32√2π²) E2/E4\n", N, L, a, h);
    const double norma = 32*sqrt(2.0)*M_PI*M_PI, tau = 0.1*h*h/(1 + 8.0/(h*h)*0.0 + 4.0);
    double E2 = 0, E4 = 0;
    for (int it = 0; it <= iters; it++) {
        // derivadas centrales
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k);
            d[0][c] = (0.5/h)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]);
            d[1][c] = (0.5/h)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]);
            d[2][c] = (0.5/h)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
        }
        // flujo de "corriente" W_i = Σ_j F_ij (∂_j n × n), y energía
        vector<V> W[3]; for (auto& v : W) v.assign(N*N*N, {0,0,0});
        E2 = 0; E4 = 0;
        #pragma omp parallel for collapse(2) reduction(+:E2,E4)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); V nn = n[c];
            double F[3][3];
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) F[p][q] = dot(nn, cross(d[p][c], d[q][c]));
            for (int p = 0; p < 3; p++) { V s = {0,0,0};
                for (int q = 0; q < 3; q++) s = s + F[p][q]*cross(d[q][c], nn);
                W[p][c] = s; E2 += dot(d[p][c], d[p][c]); }
            E4 += F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2];   // ½ Σ_ij F_ij² = Σ_{i<j} F_ij²
        }
        E2 *= h*h*h; E4 *= h*h*h;
        if (it % 200 == 0) { fprintf(f, "%d %.5f %.5f %.5f %.4f\n", it, E2, E4, (E2+E4)/norma, E2/E4); fflush(f);
            if (it % 2000 == 0) printf("it %6d  E2 %.4f  E4 %.4f  E/(32√2π²) %.5f  E2/E4 %.4f\n", it, E2, E4, (E2+E4)/norma, E2/E4); }
        if (it == iters) break;
        if (it > 0 && it % 2000 == 0) guarda();
        // gradiente: G = −2∇²n − 2 ∂_i W_i ; paso proyectado
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k);
            V lap = (1/(h*h))*(n[I(i+1,j,k)] + n[I(i-1,j,k)] + n[I(i,j+1,k)] + n[I(i,j-1,k)] + n[I(i,j,k+1)] + n[I(i,j,k-1)] - 6.0*n[c]);
            V divW = (0.5/h)*(W[0][I(i+1,j,k)] - W[0][I(i-1,j,k)] + W[1][I(i,j+1,k)] - W[1][I(i,j-1,k)] + W[2][I(i,j,k+1)] - W[2][I(i,j,k-1)]);
            V g = -2.0*lap - 2.0*divW; g = g - dot(g, n[c])*n[c];
            G[c] = g;
        }
        double dtau = (getenv("DTAU") ? atof(getenv("DTAU")) : 0.02)*h*h;
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); V m = n[c] - dtau*G[c]; n[c] = (1/sqrt(dot(m,m)))*m; }
    }
    fclose(f);
    guarda();
    printf("final: E2 %.4f  E4 %.4f  E/(32√2π²) %.5f  E2/E4 %.4f\n", E2, E4, (E2+E4)/norma, E2/E4);
}
