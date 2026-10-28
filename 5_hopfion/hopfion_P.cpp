// hopfion_P.cpp — hopfión que GIRA en el tiempo con momento angular fijo  († Claude, 28-10-2026)
//
// Pregunta: ¿puede el giro del director (y no el término de cuarto orden κ₄) impedir que el nudo se encoja?
// Precedente: Volovik y Mineev, Sov. Phys. JETP 46, 401 (1977): en ³He-A el solitón se estabiliza por su momento.
//
// Modelo [H]: además de la torsión estática en η (n = R_ẑ(kη)·n₃, masa de Kaluza–Klein, como en hopfion_eta.cpp),
//   el director gira en el tiempo alrededor de ẑ con frecuencia ω:   n(x, t) = R_ẑ(ωt)·n₃(x).
//   Ese giro tiene energía cinética T = ω²·I y momento angular S = ∂T/∂ω = 2ω·I, con  I = ∫(1 − n_z²) d³x.
//   A S fijo (S se conserva; sería el espín), la energía del giro es  T = S²/(4I),  y la frecuencia  ω = S/(2I).
//   Las configuraciones estacionarias que giran son los mínimos de (exacto para un giro alrededor de ẑ) [D]:
//
//     E_S = ∫ [ ∂n·∂n  +  C4·½F_ijF_ij  +  k²(1 − n_z²)  +  C4·k²|∇n_z|² ]  +  S²/(4I)
//            E2         E4                 Et2               Et4               Erot
//
//   C4 = coeficiente del término de cuarto orden (C4 = 1 reproduce hopfion_eta.cpp; C4 = 0 lo quita: la prueba).
//   Escalado de Derrick (tamaño λ): E2 ∝ λ, E4 ∝ 1/λ, Et2 ∝ λ³, Et4 ∝ λ, Erot ∝ 1/λ³.
//   Virial en la solución:  E2 − E4 + 3·Et2 + Et4 − 3·Erot = 0  (última columna, dividida por E).
//   Erot ∝ 1/λ³ crece al encoger: frena el colapso aunque C4 = 0 (si el resto lo permite: la simulación lo dice).
//   Condición de cola localizada (sin radiar): ω < k. Lejos del nudo, −∇²w + (k² − ω²)w = 0: con ω ≥ k la cola
//   oscila y la energía se escapa. El programa avisa si ω ≥ k. (Es la condición de los «Q-balls».)
//
// Unidades: las de hopfion_eta.cpp (longitudes en ℓ, energías en √(κ₂κ₄)) cuando C4 = 1. Con C4 = 0 no hay ℓ:
//   el tamaño lo fijan k y S, y las unidades de longitud son arbitrarias (la caja y k las fijan).
//
// Uso:  K=1 S=200 C4=0 OMP_NUM_THREADS=4 ./hopfion_P N L a iter [continuar.bin]
//   N, L: malla N³ en [−L/2, L/2]³;  a: tamaño del mapa de Hopf inicial;  iter: iteraciones.
//   Con un .bin (de hopfion.cpp, hopfion_eta.cpp o este programa, misma N y L) continúa desde él.
//   Guarda n_P_N<N>_K<k>_S<S>_C<C4>.bin cada 2000 iteraciones y al final; hopfP_...dat con la historia.
//   DTAU=... fija el paso inicial (en unidades de h²); se reduce solo si la energía sube.
// Compilar: g++ -O3 -march=native -fopenmp hopfion_P.cpp -o hopfion_P
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
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
static double env(const char* s, double def) { const char* v = getenv(s); return v ? atof(v) : def; }

int main(int argc, char** argv) {
    if (argc < 5) { fprintf(stderr, "uso: K=1 S=200 C4=0 ./hopfion_P N L a iter [continuar.bin]\n"); return 1; }
    const double KK = env("K", 1.0), k2 = KK*KK, SS = env("S", 0.0), C4 = env("C4", 1.0);
    N = atoi(argv[1]); double L = atof(argv[2]), a = atof(argv[3]); int iters = atoi(argv[4]);
    h = L/(N-1);
    vector<V> n(N*N*N), G(N*N*N); vector<V> d[3]; for (auto& v : d) v.resize(N*N*N);
    vector<V> W[3]; for (auto& v : W) v.assign(N*N*N, {0,0,0});
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        double x = (-L/2 + i*h)/a, y = (-L/2 + j*h)/a, z = (-L/2 + k*h)/a, r2 = x*x + y*y + z*z;
        double u1 = 2*x/(1+r2), u2 = 2*y/(1+r2), u3 = 2*z/(1+r2), u4 = (r2-1)/(1+r2);
        double z1r = u1, z1i = u2, z0r = u4, z0i = u3;
        double re = z1r*z0r + z1i*z0i, im = z1i*z0r - z1r*z0i, m1 = z1r*z1r + z1i*z1i, m0 = z0r*z0r + z0i*z0i;
        n[I(i,j,k)] = {2*re/(m0+m1), 2*im/(m0+m1), (m0-m1)/(m0+m1)};
        if (i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1) n[I(i,j,k)] = {0, 0, 1};
    }
    if (argc > 5) { FILE* r = fopen(argv[5], "rb"); int N2; double L2;
        if (!r) { fprintf(stderr, "no encuentro el fichero de continuación %s\n", argv[5]); return 1; }
        if (fread(&N2, sizeof N2, 1, r) != 1 || fread(&L2, sizeof L2, 1, r) != 1 || N2 != N || fabs(L2 - L) > 1e-9) {
            fprintf(stderr, "el fichero de continuación no coincide en N o L\n"); return 1; }
        if (fread(n.data(), sizeof(V), n.size(), r) != n.size()) { fprintf(stderr, "fichero incompleto\n"); return 1; }
        fclose(r); printf("continuando desde %s\n", argv[5]); }
    char tag[96]; snprintf(tag, sizeof tag, "N%d_K%g_S%g_C%g", N, KK, SS, C4);
    auto guarda = [&]() { char nb[128]; snprintf(nb, sizeof nb, "n_P_%s.bin", tag); FILE* b = fopen(nb, "wb");
        fwrite(&N, sizeof N, 1, b); fwrite(&L, sizeof L, 1, b); fwrite(n.data(), sizeof(V), n.size(), b); fclose(b); };
    char nom[128]; snprintf(nom, sizeof nom, "hopfP_%s.dat", tag); FILE* f = fopen(nom, argc > 5 ? "a" : "w");
    fprintf(f, "# N=%d L=%g a=%g h=%g K=%g S=%g C4=%g\n# iter E2 E4 Et2 Et4 Erot Etotal I omega R_rms virial\n", N, L, a, h, KK, SS, C4);
    double fmax = 1.0/(12*(1 + k2*C4 + 1e-12)), fdt = env("DTAU", std::min(0.02, 0.6*fmax));
    printf("K=%g S=%g C4=%g   paso inicial dtau = %.4f·h²\n", KK, SS, C4, fdt);
    double Eprev = 1e300; vector<V> nprev = n;
    double E2=0, E4=0, Et2=0, Et4=0, Erot=0, In=0, R2w=0, om=0;
    for (int it = 0; it <= iters; it++) {
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k);
            d[0][c] = (0.5/h)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]);
            d[1][c] = (0.5/h)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]);
            d[2][c] = (0.5/h)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
        }
        double e2=0, e4=0, it2=0, it4=0, ii=0, r2w=0;
        #pragma omp parallel for collapse(2) reduction(+:e2,e4,it2,it4,ii,r2w)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); V nn = n[c]; double F[3][3];
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) F[p][q] = dot(nn, cross(d[p][c], d[q][c]));
            for (int p = 0; p < 3; p++) { V s = {0,0,0};
                for (int q = 0; q < 3; q++) s = s + F[p][q]*cross(d[q][c], nn);
                W[p][c] = s; e2 += dot(d[p][c], d[p][c]); }
            e4 += F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2];
            double w = 1 - nn.z*nn.z; ii += w;
            it4 += d[0][c].z*d[0][c].z + d[1][c].z*d[1][c].z + d[2][c].z*d[2][c].z;
            double x = -L/2 + i*h, y = -L/2 + j*h, z = -L/2 + k*h; r2w += (x*x + y*y + z*z)*w;
        }
        double h3 = h*h*h;
        E2 = e2*h3; E4 = C4*e4*h3; In = ii*h3; Et2 = k2*In; Et4 = C4*k2*it4*h3; R2w = r2w*h3;
        if (In < 1e-12) { fprintf(stderr, "it %d: I = 0, el nudo ha desaparecido\n", it); break; }
        Erot = SS*SS/(4*In); om = SS/(2*In);
        double Et = E2 + E4 + Et2 + Et4 + Erot, vir = (E2 - E4 + 3*Et2 + Et4 - 3*Erot)/Et;
        if (it > 0 && Et > Eprev*(1 + 1e-9)) {
            fdt *= 0.5; n = nprev;
            fprintf(stderr, "it %d: la energía sube (%.4f -> %.4f); deshago el paso y bajo dtau a %.5f·h²\n", it, Eprev, Et, fdt);
            if (fdt < 1e-6) { fprintf(stderr, "dtau demasiado pequeño, paro\n"); return 2; }
            it--; continue; }
        Eprev = Et; nprev = n;
        double Rrms = sqrt(R2w/In);
        if (it % 200 == 0) { fprintf(f, "%d %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.4f %.5f\n", it, E2, E4, Et2, Et4, Erot, Et, In, om, Rrms, vir); fflush(f); }
        if (it % 2000 == 0) {
            printf("it %6d  E2 %.2f E4 %.2f Et2 %.2f Et4 %.2f Erot %.2f  E %.3f  I %.3f  ω %.4f  R_rms %.3f  virial %.4f%s\n",
                   it, E2, E4, Et2, Et4, Erot, Et, In, om, Rrms, vir, (om >= KK && SS > 0) ? "  ¡ω ≥ k: la cola radiaría!" : "");
            fflush(stdout); }
        if (it == iters) break;
        if (it > 0 && it % 2000 == 0) guarda();
        // gradiente de E_S respecto a n:  −2∇²n − 2C4·∂_iW_i + [−2k²n_z − 2C4k²∇²n_z + (S²/(2I²))·n_z] ẑ
        const double crot = SS*SS/(2*In*In);
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k);
            V lap = (1/(h*h))*(n[I(i+1,j,k)] + n[I(i-1,j,k)] + n[I(i,j+1,k)] + n[I(i,j-1,k)] + n[I(i,j,k+1)] + n[I(i,j,k-1)] - 6.0*n[c]);
            V divW = (0.5/h)*(W[0][I(i+1,j,k)] - W[0][I(i-1,j,k)] + W[1][I(i,j+1,k)] - W[1][I(i,j-1,k)] + W[2][I(i,j,k+1)] - W[2][I(i,j,k-1)]);
            double lapz = (n[I(i+1,j,k)].z + n[I(i-1,j,k)].z + n[I(i,j+1,k)].z + n[I(i,j-1,k)].z + n[I(i,j,k+1)].z + n[I(i,j,k-1)].z - 6*n[c].z)/(h*h);
            V g = -2.0*lap - (2.0*C4)*divW; g.z += -2*k2*n[c].z - 2*C4*k2*lapz + crot*n[c].z; g = g - dot(g, n[c])*n[c];
            G[c] = g;
        }
        double dtau = fdt*h*h;
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); V m = n[c] - dtau*G[c]; n[c] = (1/sqrt(dot(m,m)))*m; }
    }
    fclose(f); guarda();
    printf("final %s: E2 %.3f E4 %.3f Et2 %.3f Et4 %.3f Erot %.3f  E %.3f  I %.3f  ω %.4f (k = %g)%s\n", tag, E2, E4, Et2, Et4, Erot,
           E2+E4+Et2+Et4+Erot, In, om, KK, (om >= KK && SS > 0) ? "  ¡ω ≥ k!" : "");
    printf("Carga de Hopf: python3 hopfq.py n_P_%s.bin\n", tag);
}
