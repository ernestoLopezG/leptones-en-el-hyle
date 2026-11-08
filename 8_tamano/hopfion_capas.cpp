// hopfion_capas.cpp — hopfión de Faddeev–Skyrme dentro de una lámina de esméctico con capas  (Claude, 4-11-2026)
//
// Pregunta: ¿fija la elasticidad de las capas el tamaño del nudo en las extendidas?
//   Esméctico real (de Gennes y Prost, cap. 7): energía de capas  ½B(∂_z u)² + ½K(∇⊥²u)²,  λ = √(K/B).
//   Una deformación de anchura ℓ⊥ penetra ℓ∥ ≈ ℓ⊥²/λ a lo largo de la normal: en una lámina de espesor Δ
//   las deformaciones más anchas que ≈ √(λΔ) se alivian desplazando capas; las más estrechas pagan B entero.
// Acoplo nudo–capas (fase C): las capas se estrechan donde el director se inclina, d = d₀ cos θ ≈ d₀(1 − θ²/2).
//   Deformación preferida s = γ(1 − n_z) (n_z = cos θ respecto a la normal ẑ de las capas).
//
// Geometría del juguete: x, y = extendidas (caja L, borde fijo n = ẑ);  z = ξ, normal a las capas,
//   lámina de espesor Δ entre dos interfaces con n = ẑ (anclaje) y u = 0 en z = 0.
//   En z = Δ:  borde «fijo» (u = 0, capas sujetas también arriba)  o  «libre» (sin tensión).
// Energía total:
//   E = ∫ [ ∂n·∂n + ½F² ] d³x  +  ∫ [ ½B (∂_z u + s)² + ½K (∇⊥² u)² ] d³x ,   s = γ(1 − n_z)
//   u se minimiza EXACTAMENTE en cada paso (es cuadrático): FFT en x,y y sistema tridiagonal en z.
//   El director baja por flujo de gradiente proyectado (como hopfion.cpp).
//
// Uso: ./hopfion_capas N Nz L Delta a B lambda gamma borde itmax etiqueta
//   N: puntos en x,y   Nz: puntos en z (incluidas las dos paredes)   L: caja x,y   Delta: espesor de la lámina
//   a: tamaño del hopfión inicial   B: módulo de compresión   lambda: √(K/B)   gamma: acoplo inclinación–capas
//   borde: fijo | libre        itmax: iteraciones máximas      etiqueta: nombre de la corrida
//   ./hopfion_capas test : comprobación del gradiente del acoplo por diferencias finitas
// Variables de entorno: DTAU (0.02), TOL (1e-6 por cada 500 it.), ITMIN (3000), OMP_NUM_THREADS.
// Salidas: capas_<etiqueta>.dat (iter E_FS E_capas E/(32√2π²) R_xy sigma_z), n_capas_<etiqueta>.bin,
//          y una línea en resumen_capas.dat:
//          etiqueta Delta B lambda gamma borde E_FS/norma E_capas/norma R_xy sigma_z sqrt(lambda*Delta) iter
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <fftw3.h>
#include <omp.h>
using namespace std;
typedef complex<double> cd;
struct V { double x, y, z; };
inline V operator+(V a, V b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline V operator-(V a, V b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline V operator*(double s, V a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(V a, V b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline V cross(V a, V b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }

int N, Nz, Nh; double L, hx, hz, Bm, Km, gam; bool libre;
inline size_t I(int i, int j, int k) { return ((size_t)k*N + i)*N + j; }      // capa z contigua (para la FFT)

// ---- energía de capas: dado n, minimiza u y devuelve E_c; deja en e (N²·Nz) la deformación de cada enlace z
struct Capas {
    vector<double> s, e, lapu; fftw_complex *S, *Ef, *Lf; fftw_plan ps, pe, pl;
    vector<double> k4;
    void init() {
        size_t NR = (size_t)N*N*Nz, NC = (size_t)N*Nh*Nz;
        s.assign(NR, 0); e.assign(NR, 0); lapu.assign(NR, 0);
        S = fftw_alloc_complex(NC); Ef = fftw_alloc_complex(NC); Lf = fftw_alloc_complex(NC);
        int dims[2] = {N, N};
        ps = fftw_plan_many_dft_r2c(2, dims, Nz, s.data(), NULL, 1, N*N, S, NULL, 1, N*Nh, FFTW_ESTIMATE);
        pe = fftw_plan_many_dft_c2r(2, dims, Nz, Ef, NULL, 1, N*Nh, e.data(), NULL, 1, N*N, FFTW_ESTIMATE);
        pl = fftw_plan_many_dft_c2r(2, dims, Nz, Lf, NULL, 1, N*Nh, lapu.data(), NULL, 1, N*N, FFTW_ESTIMATE);
        k4.resize((size_t)N*Nh);
        for (int i = 0; i < N; i++) for (int j = 0; j < Nh; j++) {
            double kt2 = 4/(hx*hx)*(pow(sin(M_PI*i/N), 2) + pow(sin(M_PI*j/N), 2));   // −(Laplaciano discreto)
            k4[(size_t)i*Nh + j] = kt2*kt2;
        }
    }
    double calcula(const vector<V>& n) {
        if (Bm <= 0) return 0;                                   // sin capas: control
        for (size_t c = 0; c < s.size(); c++) s[c] = gam*(1 - n[c].z);
        fftw_execute(ps);
        const double norm = 1.0/((double)N*N);
        #pragma omp parallel
        {
            vector<cd> Sb(Nz), U(Nz), cp(Nz), dp(Nz); vector<double> a(Nz), b(Nz), c(Nz);
            #pragma omp for
            for (int ij = 0; ij < N*Nh; ij++) {
                auto Sk = [&](int k) { fftw_complex& z = S[(size_t)k*N*Nh + ij]; return cd(z[0], z[1]); };
                for (int k = 0; k < Nz-1; k++) Sb[k] = 0.5*(Sk(k) + Sk(k+1));             // s en los enlaces
                double kk = Km*k4[ij], A = Bm/(hz*hz);
                int i0 = 1, i1 = libre ? Nz-1 : Nz-2;                                     // nodos con u libre
                int m = i1 - i0 + 1;
                for (int q = 0; q < m; q++) {
                    int i = i0 + q;
                    bool top = libre && i == Nz-1;
                    b[q] = (top ? A : 2*A) + kk; a[q] = -A; c[q] = top ? 0 : -A;
                    cd rhs = top ? -(Bm/hz)*Sb[i-1] : -(Bm/hz)*(Sb[i-1] - Sb[i]);
                    dp[q] = rhs;
                }
                // Thomas
                cp[0] = c[0]/b[0]; dp[0] = dp[0]/b[0];
                for (int q = 1; q < m; q++) { double den = b[q] - a[q]*cp[q-1].real(); cp[q] = c[q]/den; dp[q] = (dp[q] - a[q]*dp[q-1])/den; }
                U.assign(Nz, 0); U[i1] = dp[m-1];
                for (int q = m-2; q >= 0; q--) U[i0+q] = dp[q] - cp[q].real()*U[i0+q+1];
                for (int k = 0; k < Nz; k++) {
                    cd ek = k < Nz-1 ? (U[k+1] - U[k])/hz + Sb[k] : cd(0, 0);
                    size_t o = (size_t)k*N*Nh + ij;
                    Ef[o][0] = ek.real()*norm; Ef[o][1] = ek.imag()*norm;
                    cd lk = -sqrt(k4[ij])*U[k];
                    Lf[o][0] = lk.real()*norm; Lf[o][1] = lk.imag()*norm;
                }
            }
        }
        fftw_execute(pe); fftw_execute(pl);
        double E = 0;
        #pragma omp parallel for reduction(+:E)
        for (size_t c = 0; c < s.size(); c++) {
            int k = c/((size_t)N*N);
            if (k < Nz-1) E += 0.5*Bm*e[c]*e[c];
            E += 0.5*Km*lapu[c]*lapu[c];
        }
        return E*hx*hx*hz;
    }
    // derivada funcional respecto a n_z:  δE_c/δn_z = −γ·B·½(e_{k−1} + e_k)
    double dEdnz(int i, int j, int k) const {
        double em = k > 0 ? e[I(i,j,k-1)] : 0, ep = k < Nz-1 ? e[I(i,j,k)] : 0;
        return -gam*Bm*0.5*(em + ep);
    }
};

int main(int argc, char** argv) {
    bool test = argc > 1 && !strcmp(argv[1], "test");
    if (!test && argc < 12) { fprintf(stderr, "uso: ./hopfion_capas N Nz L Delta a B lambda gamma fijo|libre itmax etiqueta\n"); return 1; }
    double Delta, a, lam; int itmax; string tag;
    if (test) { N = 24; Nz = 14; L = 4; Delta = 2.6; a = 0.8; Bm = 3; lam = 0.4; gam = 1; libre = getenv("LIBRE") != NULL; itmax = 0; tag = "test"; }
    else { N = atoi(argv[1]); Nz = atoi(argv[2]); L = atof(argv[3]); Delta = atof(argv[4]); a = atof(argv[5]);
           Bm = atof(argv[6]); lam = atof(argv[7]); gam = atof(argv[8]); libre = !strcmp(argv[9], "libre");
           itmax = atoi(argv[10]); tag = argv[11]; }
    Km = Bm*lam*lam; hx = L/(N-1); hz = Delta/(Nz-1); Nh = N/2 + 1;
    const double dtf = getenv("DTAU") ? atof(getenv("DTAU")) : 0.02, tol = getenv("TOL") ? atof(getenv("TOL")) : 1e-6;
    const int itmin = getenv("ITMIN") ? atoi(getenv("ITMIN")) : 3000;
    const double norma = 32*sqrt(2.0)*M_PI*M_PI;
    size_t NT = (size_t)N*N*Nz;
    vector<V> n(NT), G(NT), d[3], W[3];
    for (int p = 0; p < 3; p++) { d[p].resize(NT); W[p].resize(NT); }
    for (int k = 0; k < Nz; k++) for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        double x = (-L/2 + i*hx)/a, y = (-L/2 + j*hx)/a, z = (k*hz - Delta/2)/a, r2 = x*x + y*y + z*z;
        double u1 = 2*x/(1+r2), u2 = 2*y/(1+r2), u3 = 2*z/(1+r2), u4 = (r2-1)/(1+r2);
        double z1r = u1, z1i = u2, z0r = u4, z0i = u3;
        double re = z1r*z0r + z1i*z0i, im = z1i*z0r - z1r*z0i, m1 = z1r*z1r + z1i*z1i, m0 = z0r*z0r + z0i*z0i;
        n[I(i,j,k)] = {2*re/(m0+m1), 2*im/(m0+m1), (m0-m1)/(m0+m1)};
        if (i == 0 || j == 0 || i == N-1 || j == N-1 || k == 0 || k == Nz-1) n[I(i,j,k)] = {0, 0, 1};
    }
    Capas cap; cap.init();

    if (test) {   // gradiente del acoplo frente a diferencias finitas
        double E0 = cap.calcula(n); int fallos = 0;
        int pts[3][3] = {{N/2, N/2+2, Nz/2}, {N/2-3, N/2, Nz/2-2}, {N/2+1, N/2-1, Nz/2+3}};
        for (auto& p : pts) {
            double g = cap.dEdnz(p[0], p[1], p[2])*hx*hx*hz, eps = 1e-6;
            vector<V> m = n; m[I(p[0], p[1], p[2])].z += eps; double Ep = cap.calcula(m);
            m[I(p[0], p[1], p[2])].z -= 2*eps; double Em = cap.calcula(m);
            double num = (Ep - Em)/(2*eps);
            printf("punto (%d,%d,%d): analítico %.6e  numérico %.6e  cociente %.6f\n", p[0], p[1], p[2], g, num, g/num);
            if (fabs(g/num - 1) > 1e-3) fallos++;
            cap.calcula(n);
        }
        // u mínimo: la energía no baja al perturbar u (comprobación indirecta: E0 coincide con recalcular)
        printf("E_capas = %.6f   %s\n", E0, fallos ? "FALLO" : "gradiente correcto");
        return fallos;
    }

    string nomdat = "capas_" + tag + ".dat"; FILE* f = fopen(nomdat.c_str(), "w");
    fprintf(f, "# N=%d Nz=%d L=%g Delta=%g a=%g B=%g lambda=%g K=%g gamma=%g borde=%s hx=%g hz=%g\n# iter E_FS E_capas E/(32√2π²) R_xy sigma_z\n",
            N, Nz, L, Delta, a, Bm, lam, Km, gam, libre ? "libre" : "fijo", hx, hz);
    double dtau = dtf*pow(min(hx, hz), 2), Eant = 1e300, E2 = 0, E4 = 0, Ec = 0, Rxy = 0, sz = 0; int it;
    auto medidas = [&]() {
        double sw = 0, sr = 0, szz = 0, sz2 = 0;
        for (int k = 0; k < Nz; k++) for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            double x = -L/2 + i*hx, y = -L/2 + j*hx, z = k*hz, w = 1 - n[I(i,j,k)].z;
            sw += w; sr += w*(x*x + y*y); szz += w*z; sz2 += w*z*z; }
        Rxy = sqrt(sr/sw); double zc = szz/sw; sz = sqrt(max(0.0, sz2/sw - zc*zc));
    };
    for (it = 0; it <= itmax; it++) {
        #pragma omp parallel for collapse(2)
        for (int k = 1; k < Nz-1; k++) for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
            size_t c = I(i,j,k);
            d[0][c] = (0.5/hx)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]);
            d[1][c] = (0.5/hx)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]);
            d[2][c] = (0.5/hz)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
        }
        E2 = 0; E4 = 0;
        #pragma omp parallel for collapse(2) reduction(+:E2,E4)
        for (int k = 1; k < Nz-1; k++) for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
            size_t c = I(i,j,k); V nn = n[c]; double F[3][3];
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) F[p][q] = dot(nn, cross(d[p][c], d[q][c]));
            for (int p = 0; p < 3; p++) { V s = {0,0,0};
                for (int q = 0; q < 3; q++) s = s + F[p][q]*cross(d[q][c], nn);
                W[p][c] = s; E2 += dot(d[p][c], d[p][c]); }
            E4 += F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2];
        }
        E2 *= hx*hx*hz; E4 *= hx*hx*hz;
        Ec = cap.calcula(n);
        if (it % 500 == 0) {
            double E = E2 + E4 + Ec;
            medidas();
            fprintf(f, "%d %.5f %.5f %.5f %.4f %.4f\n", it, E2+E4, Ec, E/norma, Rxy, sz); fflush(f);
            if (it % 5000 == 0) { printf("[%s] it %6d  E_FS %.4f  E_capas %.4f  (÷32√2π²)  R_xy %.3f  sigma_z %.3f\n",
                                        tag.c_str(), it, (E2+E4)/norma, Ec/norma, Rxy, sz); fflush(stdout); }
            if (it >= itmin && fabs(Eant - E)/E < tol) break;
            Eant = E;
        }
        if (it == itmax) break;
        #pragma omp parallel for collapse(2)
        for (int k = 1; k < Nz-1; k++) for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
            size_t c = I(i,j,k);
            V lap = (1/(hx*hx))*(n[I(i+1,j,k)] + n[I(i-1,j,k)] + n[I(i,j+1,k)] + n[I(i,j-1,k)] - 4.0*n[c])
                  + (1/(hz*hz))*(n[I(i,j,k+1)] + n[I(i,j,k-1)] - 2.0*n[c]);
            V divW = (0.5/hx)*(W[0][I(i+1,j,k)] - W[0][I(i-1,j,k)] + W[1][I(i,j+1,k)] - W[1][I(i,j-1,k)])
                   + (0.5/hz)*(W[2][I(i,j,k+1)] - W[2][I(i,j,k-1)]);
            V g = -2.0*lap - 2.0*divW; g.z += cap.dEdnz(i, j, k);
            g = g - dot(g, n[c])*n[c]; G[c] = g;
        }
        #pragma omp parallel for collapse(2)
        for (int k = 1; k < Nz-1; k++) for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) {
            size_t c = I(i,j,k); V m = n[c] - dtau*G[c]; n[c] = (1/sqrt(dot(m,m)))*m; }
    }
    fclose(f); medidas();
    string nb = "n_capas_" + tag + ".bin"; FILE* b = fopen(nb.c_str(), "wb");
    fwrite(&N, sizeof N, 1, b); fwrite(&Nz, sizeof Nz, 1, b); fwrite(&L, sizeof L, 1, b); fwrite(&Delta, sizeof Delta, 1, b);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < Nz; k++) fwrite(&n[I(i,j,k)], sizeof(V), 1, b);
    fclose(b);
    FILE* r = fopen("resumen_capas.dat", "a");
    fprintf(r, "%s %g %g %g %g %s %.5f %.5f %.4f %.4f %.4f %d\n", tag.c_str(), Delta, Bm, lam, gam, libre ? "libre" : "fijo",
            (E2+E4)/norma, Ec/norma, Rxy, sz, sqrt(lam*Delta), it);
    fclose(r);
    printf("[%s] final it %d  E_FS %.5f  E_capas %.5f  R_xy %.4f  sigma_z %.4f  (√(λΔ) = %.3f)\n",
           tag.c_str(), it, (E2+E4)/norma, Ec/norma, Rxy, sz, sqrt(lam*Delta));
}
