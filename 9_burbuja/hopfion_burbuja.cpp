// hopfion_burbuja.cpp — hopfión dentro de una burbuja ordenada en un fondo desordenado  (Claude, 10-11-2026)
//
// Pregunta: si el director de las extendidas no tiene una dirección global (fondo desordenado) y el orden
//   solo existe en una burbuja alrededor del nudo (inducida por la pulsación, sección 4), ¿sobrevive el nudo?
//   ¿Qué radio mínimo de burbuja necesita, comparado con su propio tamaño?
//
// Modelo: grado de orden S(r) impuesto (no evoluciona), S = 1 dentro, cae a 0 fuera:
//     S(r) = ½[1 − tanh((r − Rb)/w)]
//   Energía de Faddeev–Skyrme pesada por el orden (la rigidez de un cristal líquido va como S², y el término
//   de cuarto orden, de cuatro derivadas, como S⁴):
//     E = ∫ [ S² ∂_i n·∂_i n  +  ½ S⁴ F_ij F_ij ] d³x ,   F_ij = n·(∂_i n × ∂_j n)
//   Donde S = 0 el director no cuesta nada: es el fondo desordenado. No hay borde fijo.
//   Flujo de gradiente proyectado:  G = −2 ∂_i(S² ∂_i n) − 2 ∂_i(S⁴ W_i),   W_i = Σ_j F_ij (∂_j n × n).
// Arranque: mapa de Hopf de carga 1 y tamaño a dentro de la burbuja.  Fuera de la burbuja (r > Rb + 2w):
//     fondo = uniforme  → n = ẑ (el mismo valor que el mapa de Hopf en el infinito)
//     fondo = aleatorio → direcciones al azar en cada nodo (fondo sin orden), unidas suavemente al nudo
// Diagnóstico principal: la energía. Si el nudo se escapa hacia la zona sin orden, se ensancha y su energía
//   cae hacia 0 (Derrick). Si sobrevive, la energía se queda cerca de la del fondo ordenado (≈ 1,22–1,23).
//   Q con ventana: hopfq_burbuja.py (la Q completa no sirve con fondo aleatorio).
//
// Uso: ./hopfion_burbuja N L a Rb w fondo itmax etiqueta [semilla]
//   N: puntos por lado   L: caja [−L/2, L/2]³   a: tamaño del mapa de Hopf inicial   Rb: radio de la burbuja
//   w: anchura del borde   fondo: uniforme | aleatorio   itmax: iteraciones   etiqueta: nombre de la corrida
// Variables de entorno: DTAU (0.02), OMP_NUM_THREADS.
// Salidas: burbuja_<etiqueta>.dat (iter E2 E4 E/(32√2π²) E2/E4 fracción_dentro r_rms),
//          n_<etiqueta>.bin (N, L, Rb, w y el campo) y una línea en resumen_burbuja.dat:
//          etiqueta N L a Rb w fondo E/(32√2π²) E2/E4 fracción_dentro r_rms iter
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <random>
#include <string>
#include <omp.h>
using namespace std;
struct V { double x, y, z; };
inline V operator+(V a, V b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline V operator-(V a, V b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline V operator*(double s, V a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(V a, V b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline V cross(V a, V b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
inline V unit(V a) { return (1/sqrt(dot(a,a)))*a; }
int N; double h;
inline size_t I(int i, int j, int k) { return ((size_t)i*N + j)*N + k; }

int main(int argc, char** argv) {
    if (argc < 9) { fprintf(stderr, "uso: ./hopfion_burbuja N L a Rb w uniforme|aleatorio itmax etiqueta [semilla]\n"); return 1; }
    N = atoi(argv[1]); double L = atof(argv[2]), a = atof(argv[3]), Rb = atof(argv[4]), w = atof(argv[5]);
    bool aleat = !strcmp(argv[6], "aleatorio"); int itmax = atoi(argv[7]); string et = argv[8];
    unsigned semilla = argc > 9 ? atoi(argv[9]) : 12345;
    h = L/(N-1); size_t NT = (size_t)N*N*N;
    vector<V> n(NT), G(NT), d[3], W[3]; for (int p = 0; p < 3; p++) { d[p].resize(NT); W[p].resize(NT); }
    vector<double> S2(NT), S4(NT);
    mt19937 rng(semilla); normal_distribution<double> gau(0, 1);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        double X = -L/2 + i*h, Y = -L/2 + j*h, Z = -L/2 + k*h, r = sqrt(X*X + Y*Y + Z*Z);
        double S = 0.5*(1 - tanh((r - Rb)/w)); S2[I(i,j,k)] = S*S; S4[I(i,j,k)] = S*S*S*S;
        double x = X/a, y = Y/a, z = Z/a, r2 = x*x + y*y + z*z;
        double u1 = 2*x/(1+r2), u2 = 2*y/(1+r2), u3 = 2*z/(1+r2), u4 = (r2-1)/(1+r2);
        double re = u1*u4 + u2*u3, im = u2*u4 - u1*u3, m1 = u1*u1 + u2*u2, m0 = u4*u4 + u3*u3;
        V hopf = {2*re/(m0+m1), 2*im/(m0+m1), (m0-m1)/(m0+m1)};
        V fuera = {0, 0, 1};
        if (aleat) fuera = unit({gau(rng), gau(rng), gau(rng)});
        double t = 0.5*(1 + tanh((r - (Rb + w))/w));           // 0 dentro, 1 en el fondo
        n[I(i,j,k)] = unit((1 - t)*hopf + t*fuera);
    }
    const double norma = 32*sqrt(2.0)*M_PI*M_PI;
    double dtau = (getenv("DTAU") ? atof(getenv("DTAU")) : 0.02)*h*h;
    string nd = "burbuja_" + et + ".dat"; FILE* f = fopen(nd.c_str(), "w");
    fprintf(f, "# N=%d L=%g a=%g Rb=%g w=%g fondo=%s h=%g\n# iter E2 E4 E/(32√2π²) E2/E4 fraccion_dentro(r<Rb) r_rms\n",
            N, L, a, Rb, w, argv[6], h);
    auto guarda = [&]() { string nb = "n_" + et + ".bin"; FILE* b = fopen(nb.c_str(), "wb");
        fwrite(&N, sizeof N, 1, b); fwrite(&L, sizeof L, 1, b); fwrite(&Rb, sizeof Rb, 1, b); fwrite(&w, sizeof w, 1, b);
        fwrite(n.data(), sizeof(V), NT, b); fclose(b); };
    double E2 = 0, E4 = 0, Edentro = 0, Er2 = 0; int it;
    for (it = 0; it <= itmax; it++) {
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k);
            d[0][c] = (0.5/h)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]);
            d[1][c] = (0.5/h)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]);
            d[2][c] = (0.5/h)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
        }
        E2 = 0; E4 = 0; Edentro = 0; Er2 = 0;
        #pragma omp parallel for collapse(2) reduction(+:E2,E4,Edentro,Er2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); V nn = n[c]; double F[3][3], e2 = 0;
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) F[p][q] = dot(nn, cross(d[p][c], d[q][c]));
            for (int p = 0; p < 3; p++) { V s = {0,0,0};
                for (int q = 0; q < 3; q++) s = s + F[p][q]*cross(d[q][c], nn);
                W[p][c] = S4[c]*s; e2 += dot(d[p][c], d[p][c]); }
            double e = S2[c]*e2 + S4[c]*(F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2]);
            E2 += S2[c]*e2; E4 += e - S2[c]*e2;
            double X = -L/2 + i*h, Y = -L/2 + j*h, Z = -L/2 + k*h, r2 = X*X + Y*Y + Z*Z;
            if (r2 < Rb*Rb) Edentro += e;
            Er2 += e*r2;
        }
        E2 *= h*h*h; E4 *= h*h*h; Edentro *= h*h*h; Er2 *= h*h*h;
        double Et = E2 + E4, frac = Et > 0 ? Edentro/Et : 0, rrms = Et > 0 ? sqrt(Er2/Et) : 0;
        if (it % 200 == 0) { fprintf(f, "%d %.5f %.5f %.5f %.4f %.4f %.4f\n", it, E2, E4, Et/norma, E2/E4, frac, rrms); fflush(f); }
        if (it % 2000 == 0) printf("[%s] it %6d  E/(32√2π²) %.5f  E2/E4 %.4f  dentro %.3f  r_rms %.3f\n",
                                    et.c_str(), it, Et/norma, E2/E4, frac, rrms);
        if (it == itmax) break;
        if (it > 0 && it % 5000 == 0) guarda();
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k);
            // ∂_i(S² ∂_i n) en forma de flujo, con S² en los puntos medios
            V lap = {0,0,0};
            size_t nb[6] = {I(i+1,j,k), I(i-1,j,k), I(i,j+1,k), I(i,j-1,k), I(i,j,k+1), I(i,j,k-1)};
            for (int m = 0; m < 6; m++) lap = lap + (0.5*(S2[nb[m]] + S2[c]))*(n[nb[m]] - n[c]);
            lap = (1/(h*h))*lap;
            V divW = (0.5/h)*(W[0][I(i+1,j,k)] - W[0][I(i-1,j,k)] + W[1][I(i,j+1,k)] - W[1][I(i,j-1,k)]
                              + W[2][I(i,j,k+1)] - W[2][I(i,j,k-1)]);
            V g = -2.0*lap - 2.0*divW; g = g - dot(g, n[c])*n[c];
            G[c] = g;
        }
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); n[c] = unit(n[c] - dtau*G[c]); }
    }
    fclose(f); guarda();
    double Et = E2 + E4;
    FILE* r = fopen("resumen_burbuja.dat", "a");
    fprintf(r, "%s %d %g %g %g %g %s %.5f %.4f %.4f %.4f %d\n", et.c_str(), N, L, a, Rb, w, argv[6], Et/norma, E2/E4,
            Et > 0 ? Edentro/Et : 0, Et > 0 ? sqrt(Er2/Et) : 0, it);
    fclose(r);
    printf("[%s] final: E/(32√2π²) %.5f  E2/E4 %.4f\n", et.c_str(), Et/norma, E2/E4);
}
