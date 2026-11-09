// nucleacion_persistente.cpp — extensión de nucleacion.cpp (Claude, 7-11-2026)
//
// NOVEDADES RESPECTO A nucleacion.cpp (el resto del código es idéntico y está comprobado):
//   1) NUCLEACIÓN PERSISTENTE: además de la fusión inicial, un núcleo fundido que no se enfría nunca,
//        T_per(x) = Tp · perfil(x),   forma=bola: exp(−r²/Rp²)   forma=toro: exp(−d²/Rp²), d = distancia al anillo de radio Rt
//      Representa la nucleación que la pulsación mantiene mientras la partícula existe. Con Tp > 1 el centro del
//      perfil está fundido (S = 0) y lleva ruido térmico siempre.
//   2) k4: multiplica el término de cuarto orden (k4 = 0 lo quita). Prueba «núcleo fundido en lugar de κ₄»:
//      sin E4, el teorema de Derrick hace que el nudo se encoja hasta la malla; ¿lo impide un núcleo fundido fijo?
//   3) init=<fichero.bin>: arranca desde un director guardado (p. ej. un nudo Q = 1 ya relajado).
//   4) Medidas nuevas: R_rms = √⟨ρ²⟩ con peso 1 − n_z (tamaño), y Q_S: la misma fórmula de Q con B·S²
//      (descuenta la zona fundida, donde el ruido hace que Q no signifique nada). Ninguna de las dos es un
//      invariante topológico mientras haya zona fundida: son indicadores.
//
// --- cabecera original de nucleacion.cpp ---
// nucleacion.cpp — ¿Puede una nucleación (fusión local del orden del hyle) dejar atrapado un nudo de Hopf?
//                  (Claude, 2-11-2026; revisado y ejecutado por el autor)
//
// IDEA FÍSICA
//   La carga de Hopf Q no cambia mientras el director esté definido en todos los puntos (topología).
//   Para crear un nudo desde el vacío (Q = 0) el orden tiene que fundirse en alguna región: la pulsación,
//   al superar el umbral del sustrato, produce una nucleación (región isótropa); al volver a ordenarse,
//   el director puede quedar enlazado. Precedentes: mecanismo de Kibble–Zurek (³He, Ruutu et al. 1996;
//   Bäuerle et al. 1996) y hopfiones creados en cristales líquidos quirales con una región isótropa
//   inducida por luz (Shvetsov et al. 2025).
//   El flujo que la pulsación induce en el hyle se representa como un anillo de vórtice con remolino,
//   de helicidad de signo s. Por Mermin–Ho la vorticidad es B (el retorcimiento del director), así que el
//   flujo favorece las texturas cuya B se alinea con él: término −λ ∫ Ω·B.
//
// MODELO (límite de módulo rígido)
//   Director unitario n(x), como en hopfion.cpp. La fusión entra por el grado de orden local S(x,t):
//     S² = max(0, 1 − T(x,t)),   T = "temperatura" de la nucleación (T > 1: fundido, S = 0)
//   Es el límite β → ∞ de un parámetro de orden m = S·n con módulo rígido: (∂m)² → S²(∂n)², F(m) → S³F(n).
//   (Una primera versión con |m| libre dejaba que el nudo se deshiciera por una sola celda de la malla,
//    un atajo numérico sin física; por eso se usa este límite.)
//
// ENERGÍA (mismas unidades y convención que hopfion.cpp: longitudes en ℓ, energías en √(κ₂κ₄))
//   E = ∫ [ S² ∂_i n·∂_i n  +  S⁶ ½ F_ij F_ij  −  λ S³ Ω·B ] d³x
//   F_ij = n·(∂_i n × ∂_j n),   B = (F_yz, F_zx, F_xy)   (vorticidad de Mermin–Ho)
//   Donde S = 0 el director no tiene rigidez y el ruido térmico lo desordena: es la región fundida.
//
// DINÁMICA: flujo de gradiente proyectado sobre la esfera + ruido térmico (solo donde T > 0).
//   Fases:  (1) fusión: T = T0·exp(−r²/Rc²) durante thold
//           (2) enfriamiento: T decae como exp(−(t−thold)/tauq); el flujo (λ) sigue actuando
//           (3) retirada del flujo: λ baja linealmente a 0 en trampa; después, relajación libre durante relax
//               (λ = 0, T = 0, sin ruido)  → ¿sobrevive el nudo por sí solo?
//           Con retira=0 el flujo se mantiene hasta el final: es el caso del modelo, en el que la pulsación
//           (y su flujo) existe mientras existe la partícula; se pregunta entonces si el nudo que sostiene es estable.
//   Borde: n = ẑ fijo.
//
// MEDIDAS: Q por la fórmula de hopfq.py, Q = (1/16π²) ∫ A·B, con A en el gauge A_z = 0
//   (A_x = ∫ B_y dz, A_y = −∫ B_x dz). Durante la fusión Q no es un entero (el director está desordenado).
//   Referencia: el hopfión relajado tiene E2 + E4 ≈ 544 (sección 4 del preprint; caja 12, malla 150³);
//   en cajas pequeñas la energía es mayor porque el borde lo comprime.
//
// COMPROBACIONES (./nucleacion test):
//   1) Q del mapa de Hopf analítico ≈ −1 (comprueba la fórmula de Q; mismo signo que hopfion.cpp).
//   2) Derivada direccional: [E(n) − E(n − εG)] / (ε ∫G²) ≈ 1 para cada término (comprueba el gradiente).
//
// USO:  ./nucleacion clave=valor ...      (ver los valores por defecto abajo)
//   ej.: ./nucleacion N=96 L=8 lam=2 s=1 seed=1 tag=s+1_a
//   control de estabilidad: ./nucleacion init=hopf T0=0 lam=0 thold=0 tauq=0 trampa=0 relax=1 tag=ctl
// SALIDAS: log_<tag>.dat (serie temporal), corte_<tag>_<paso>.dat (plano y = 0: x z S n_z),
//          n_<tag>.bin (director final, mismo formato que hopfion.cpp → vale para hopfq.py y figura_hopfion.py),
//          una línea por ejecución en resumen.dat
//
// Compilar:  g++ -O3 -march=native -fopenmp nucleacion.cpp -o nucleacion
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <omp.h>
using namespace std;
struct V { double x, y, z; };
inline V operator+(V a, V b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
inline V operator-(V a, V b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
inline V operator*(double s, V a) { return {s*a.x, s*a.y, s*a.z}; }
inline double dot(V a, V b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline V cross(V a, V b) { return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
inline V unit(V a) { double r = sqrt(dot(a,a)); return r > 1e-300 ? (1/r)*a : V{0,0,1}; }

// ---------------- parámetros ----------------
int N = 96; double L = 8.0, h;
double T0 = 2.0;         // temperatura máxima de la nucleación (T > 1 = fundido)
double Rc = 1.4;         // radio de la nucleación
double thold = 0.10;     // duración de la fusión
double tauq = 0.15;      // tiempo de enfriamiento (se enfría durante 5·tauq)
double theta = 0.05;     // intensidad del ruido térmico en la zona caliente
double lam = 2;          // intensidad del flujo inducido por la pulsación (0 = sin flujo)
double s = 1;            // helicidad del flujo: +1, −1, 0 (anillo sin remolino)
double R0 = 1.2, w = 0.6;// radio y grosor del anillo de vórtice del flujo
double trampa = 0.3;     // al acabar el enfriamiento, λ baja linealmente a 0 en este tiempo
double relax = 1.0;      // duración de la relajación libre final
int retira = 1;          // 1: se retira el flujo al final; 0: el flujo se mantiene (la pulsación sigue ahí)
double dtf = 0.02;       // paso de tiempo = dtf·h²  (el mismo de hopfion.cpp)
double Tp = 0, Rp = 0.35, Rt = 1.0, k4 = 1; string forma = "toro";
int seed = 1; string tag = "prueba", init = "vacio"; double a_hopf = 1.0;
int cada = 50, cortes = 2000;

inline size_t I(int i, int j, int k) { return ((size_t)i*N + j)*N + k; }
inline double X(int i) { return -L/2 + i*h; }

// ---------------- campos auxiliares ----------------
vector<V> Om;            // vorticidad del flujo inducido Ω(x)
vector<double> S2;       // S² = max(0, 1 − T)
vector<V> d[3], J[3], G;

void prepara_omega() {
    Om.assign((size_t)N*N*N, {0,0,0});
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        double x = X(i), y = X(j), z = X(k), rho = sqrt(x*x + y*y);
        if (rho < 1e-12) continue;
        V er = {x/rho, y/rho, 0}, ephi = {-y/rho, x/rho, 0}, ez = {0, 0, 1};
        double dd = sqrt((rho-R0)*(rho-R0) + z*z), g = exp(-dd*dd/(w*w));
        V epol = dd > 1e-12 ? (1/dd)*((-z)*er + (rho-R0)*ez) : V{0,0,0};
        Om[I(i,j,k)] = g*(ephi + s*epol);   // a lo largo del núcleo del anillo (φ̂) + remolino alrededor del núcleo (s)
    }
}
double Tper(int i, int j, int k) {           // núcleo persistente
    if (Tp <= 0) return 0;
    double x = X(i), y = X(j), z = X(k);
    if (forma == "bola") return Tp*exp(-(x*x + y*y + z*z)/(Rp*Rp));
    double rho = sqrt(x*x + y*y), d2 = (rho-Rt)*(rho-Rt) + z*z; return Tp*exp(-d2/(Rp*Rp)); }
double Tlocal(int i, int j, int k, double Tmax) {
    double r2 = X(i)*X(i) + X(j)*X(j) + X(k)*X(k); return Tmax*exp(-r2/(Rc*Rc)) + Tper(i, j, k); }
void pon_temperatura(double Tmax) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++)
        S2[I(i,j,k)] = max(0.0, 1 - Tlocal(i, j, k, Tmax));
}

struct Es { double E2, E4, EB; };
// Energías y gradiente G = δE/δn (sin proyectar; se proyecta al actualizar). solo: −1 todos; 0 E2; 1 E4; 2 EB
Es energia_y_gradiente(const vector<V>& n, double lamE, bool calculaG, int solo = -1) {
    bool t2 = solo < 0 || solo == 0, t4 = solo < 0 || solo == 1, tB = (solo < 0 || solo == 2) && lamE != 0;
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
        size_t c = I(i,j,k);
        d[0][c] = (0.5/h)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]);
        d[1][c] = (0.5/h)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]);
        d[2][c] = (0.5/h)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
    }
    double E2 = 0, E4 = 0, EB = 0;
    #pragma omp parallel for collapse(2) reduction(+:E2,E4,EB)
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
        size_t c = I(i,j,k); V nn = n[c]; double sq = S2[c], w4 = k4*sq*sq*sq, wB = sq*sqrt(sq);
        V cx[3][3];
        for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) cx[p][q] = cross(d[p][c], d[q][c]);
        double F[3][3];
        for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) F[p][q] = dot(nn, cx[p][q]);
        V nxd[3]; for (int q = 0; q < 3; q++) nxd[q] = cross(nn, d[q][c]);
        V Sl = {0,0,0}; V Jl[3] = {{0,0,0},{0,0,0},{0,0,0}};
        if (t2) for (int p = 0; p < 3; p++) E2 += sq*dot(d[p][c], d[p][c]);
        if (t4) {
            E4 += w4*(F[0][1]*F[0][1] + F[0][2]*F[0][2] + F[1][2]*F[1][2]);           // S⁶ · ½ Σ_ij F_ij²
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) Sl = Sl + (w4*F[p][q])*cx[p][q];
            for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++) Jl[p] = Jl[p] + (2*w4*F[p][q])*nxd[q];
        }
        if (tB) {
            V O = wB*Om[c]; double Bv[3] = {F[1][2], F[2][0], F[0][1]};
            EB += O.x*Bv[0] + O.y*Bv[1] + O.z*Bv[2];
            // δ(∫Ω·B)/δn = Σ_i Ω_i C_i + ∂_j Q_j ;  C_x = ∂_y n×∂_z n, C_y = ∂_z n×∂_x n, C_z = ∂_x n×∂_y n
            //   Q_j = Σ_{i,k} ε_ijk Ω_i (n × ∂_k n)        (comprobado con ./nucleacion test)
            V Cv = O.x*cx[1][2] + O.y*cx[2][0] + O.z*cx[0][1];
            Sl = Sl - lamE*Cv;
            double Oa[3] = {O.x, O.y, O.z};
            for (int jj = 0; jj < 3; jj++) { V q = {0,0,0};
                for (int ii = 0; ii < 3; ii++) for (int kk = 0; kk < 3; kk++) {
                    int e = (ii==jj||jj==kk||ii==kk) ? 0 : (((jj-ii+3)%3 == 1) ? 1 : -1);  // ε_ijk
                    if (e) q = q + (e*Oa[ii])*nxd[kk]; }
                Jl[jj] = Jl[jj] - lamE*q; }
        }
        if (calculaG) { G[c] = Sl; for (int p = 0; p < 3; p++) J[p][c] = Jl[p]; }
    }
    double h3 = h*h*h;
    if (calculaG) {
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); V g = G[c];
            if (t2) {   // −2 ∂_k(S² ∂_k n), con S² en las caras (gradiente exacto de la energía con diferencias hacia delante)
                size_t nb[6] = {I(i+1,j,k), I(i-1,j,k), I(i,j+1,k), I(i,j-1,k), I(i,j,k+1), I(i,j,k-1)};
                V lap = {0,0,0}; for (int q = 0; q < 6; q++) lap = lap + (0.5*(S2[c] + S2[nb[q]]))*(n[nb[q]] - n[c]);
                g = g - (2/(h*h))*lap; }
            auto Jv = [&](int p, int a, int b, int cc) { if (a<1||b<1||cc<1||a>N-2||b>N-2||cc>N-2) return V{0,0,0}; return J[p][I(a,b,cc)]; };
            V div = (0.5/h)*(Jv(0,i+1,j,k) - Jv(0,i-1,j,k) + Jv(1,i,j+1,k) - Jv(1,i,j-1,k) + Jv(2,i,j,k+1) - Jv(2,i,j,k-1));
            G[c] = g + div;
        }
    }
    return {E2*h3, E4*h3, EB*h3};
}

// Carga de Hopf (pesoS = true: B·S², descuenta la zona fundida)
double cargaHopf(const vector<V>& n, bool pesoS = false) {
    vector<V> B(n.size(), {0,0,0});
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
        size_t c = I(i,j,k);
        V dx = (0.5/h)*(n[I(i+1,j,k)] - n[I(i-1,j,k)]), dy = (0.5/h)*(n[I(i,j+1,k)] - n[I(i,j-1,k)]), dz = (0.5/h)*(n[I(i,j,k+1)] - n[I(i,j,k-1)]);
        B[c] = {dot(n[c], cross(dy, dz)), dot(n[c], cross(dz, dx)), dot(n[c], cross(dx, dy))};
        if (pesoS) B[c] = S2[c]*B[c];
    }
    double Q = 0;
    #pragma omp parallel for collapse(2) reduction(+:Q)
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        double Ax = 0, Ay = 0;                       // A_x = ∫ B_y dz,  A_y = −∫ B_x dz  (trapecios)
        for (int k = 0; k < N; k++) { size_t c = I(i,j,k);
            if (k > 0) { size_t p = I(i,j,k-1); Ax += 0.5*h*(B[p].y + B[c].y); Ay -= 0.5*h*(B[p].x + B[c].x); }
            Q += Ax*B[c].x + Ay*B[c].y; }
    }
    return Q*h*h*h/(16*M_PI*M_PI);
}

double radio_rms(const vector<V>& n) {
    double a = 0, b = 0;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        double w = 1 - n[I(i,j,k)].z; a += w*(X(i)*X(i) + X(j)*X(j)); b += w; }
    return b > 0 ? sqrt(a/b) : 0; }
void guarda_bin(const vector<V>& n, const string& nom) {
    FILE* b = fopen(nom.c_str(), "wb"); fwrite(&N, sizeof N, 1, b); fwrite(&L, sizeof L, 1, b);
    fwrite(n.data(), sizeof(V), n.size(), b); fclose(b);
}
void guarda_corte(const vector<V>& n, int paso, double t) {
    char nom[256]; snprintf(nom, sizeof nom, "corte_%s_%07d.dat", tag.c_str(), paso);
    FILE* f = fopen(nom, "w"); fprintf(f, "# t=%.5f  plano y=0: x z S n_z\n", t);
    int j = N/2;
    for (int i = 0; i < N; i++) for (int k = 0; k < N; k++) { size_t c = I(i,j,k);
        fprintf(f, "%.4f %.4f %.5f %.5f\n", X(i), X(k), sqrt(S2[c]), n[c].z); }
    fclose(f);
}
void reserva() { size_t M = (size_t)N*N*N; for (int p = 0; p < 3; p++) { d[p].assign(M, {0,0,0}); J[p].assign(M, {0,0,0}); }
    G.assign(M, {0,0,0}); S2.assign(M, 1.0); }

vector<V> mapa_hopf(double a) {
    vector<V> n((size_t)N*N*N);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {
        double x = X(i)/a, y = X(j)/a, z = X(k)/a, r2 = x*x + y*y + z*z;
        double u1 = 2*x/(1+r2), u2 = 2*y/(1+r2), u3 = 2*z/(1+r2), u4 = (r2-1)/(1+r2);
        double z1r = u1, z1i = u2, z0r = u4, z0i = u3;
        double re = z1r*z0r + z1i*z0i, im = z1i*z0r - z1r*z0i, m1 = z1r*z1r + z1i*z1i, m0 = z0r*z0r + z0i*z0i;
        n[I(i,j,k)] = {2*re/(m0+m1), 2*im/(m0+m1), (m0-m1)/(m0+m1)};
        if (i == 0 || j == 0 || k == 0 || i == N-1 || j == N-1 || k == N-1) n[I(i,j,k)] = {0, 0, 1};
    }
    return n;
}

int prueba() {
    printf("== Comprobación 1: Q del mapa de Hopf analítico (mismo mapa que hopfion.cpp)\n");
    for (int NN : {48, 72}) { N = NN; L = 8; h = L/(N-1); reserva();
        auto n = mapa_hopf(1.0); double Q = cargaHopf(n);
        printf("   N=%d L=%g: Q = %+.4f   (→ −1 al afinar la malla; hopfion.cpp da Q≈−0,98)\n", N, L, Q); }
    printf("== Comprobación 2: gradiente frente a derivada direccional de la energía discreta\n");
    N = 48; L = 6; h = L/(N-1); reserva(); s = 1; R0 = 1.2; w = 0.6; prepara_omega();
    auto n = mapa_hopf(1.0);
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) {   // S² suave, no uniforme
        double r2 = X(i)*X(i) + X(j)*X(j) + X(k)*X(k); S2[I(i,j,k)] = 1 - 0.6*exp(-r2/2.0); }
    auto E2f = [&](const vector<V>& m) { double s2 = 0;    // energía de gradiente con diferencias hacia delante
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) for (int k = 0; k < N; k++) { size_t c = I(i,j,k);
            if (i < N-1) { size_t a = I(i+1,j,k); V e = m[a]-m[c]; s2 += 0.5*(S2[a]+S2[c])*dot(e,e); }
            if (j < N-1) { size_t a = I(i,j+1,k); V e = m[a]-m[c]; s2 += 0.5*(S2[a]+S2[c])*dot(e,e); }
            if (k < N-1) { size_t a = I(i,j,k+1); V e = m[a]-m[c]; s2 += 0.5*(S2[a]+S2[c])*dot(e,e); } }
        return s2*h; };
    const char* nom[3] = {"E2 (gradiente, S²)", "E4 (cuarto orden, S⁶)", "EB (flujo, −λS³Ω·B)"};
    for (int t = 0; t < 3; t++) {
        double lamE = 1.0;
        Es e0 = energia_y_gradiente(n, lamE, true, t);
        double g2 = 0; vector<V> Gt = G;
        for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
            size_t c = I(i,j,k); Gt[c] = G[c] - dot(G[c], n[c])*n[c]; g2 += dot(Gt[c], Gt[c]); }   // dirección tangente
        g2 *= h*h*h;
        printf("   %-26s", nom[t]);
        for (double eps : {1e-3, 1e-4, 1e-5}) {
            vector<V> m2 = n; for (size_t c = 0; c < n.size(); c++) m2[c] = n[c] - eps*Gt[c];
            Es e1 = energia_y_gradiente(m2, lamE, false, t);
            double dE = t==0 ? E2f(n) - E2f(m2) : t==1 ? e0.E4 - e1.E4 : -lamE*(e0.EB - e1.EB);
            printf("  ε=%.0e: %.4f", eps, dE/(eps*g2)); }
        printf("   (≈ 1)\n");
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc > 1 && !strcmp(argv[1], "test")) return prueba();
    for (int a = 1; a < argc; a++) { string kv = argv[a]; size_t p = kv.find('='); if (p == string::npos) continue;
        string k = kv.substr(0, p), v = kv.substr(p+1);
        if (k=="N") N = stoi(v); else if (k=="L") L = stod(v);
        else if (k=="T0") T0 = stod(v); else if (k=="Rc") Rc = stod(v); else if (k=="thold") thold = stod(v);
        else if (k=="tauq") tauq = stod(v); else if (k=="theta") theta = stod(v); else if (k=="lam") lam = stod(v);
        else if (k=="s") s = stod(v); else if (k=="R0") R0 = stod(v); else if (k=="w") w = stod(v);
        else if (k=="trampa") trampa = stod(v); else if (k=="relax") relax = stod(v); else if (k=="dtf") dtf = stod(v);
        else if (k=="retira") retira = stoi(v);
        else if (k=="Tp") Tp = stod(v); else if (k=="Rp") Rp = stod(v); else if (k=="Rt") Rt = stod(v);
        else if (k=="k4") k4 = stod(v); else if (k=="forma") forma = v;
        else if (k=="seed") seed = stoi(v); else if (k=="init") init = v; else if (k=="a") a_hopf = stod(v);
        else if (k=="tag") tag = v; else if (k=="cada") cada = stoi(v); else if (k=="cortes") cortes = stoi(v);
        else { fprintf(stderr, "parámetro desconocido: %s\n", k.c_str()); return 1; } }
    h = L/(N-1); reserva(); prepara_omega();
    if (h > 0.09) printf("AVISO: h = %.3f > 0,09; con mallas tan gruesas el nudo se deshace solo (Apéndice A)\n", h);
    vector<V> n((size_t)N*N*N, {0,0,1});            // vacío: director uniforme, Q = 0
    if (init == "hopf") n = mapa_hopf(a_hopf);       // control: ¿sobrevive un nudo ya formado en esta malla y caja?
    else if (init.size() > 4 && init.substr(init.size()-4) == ".bin") {
        FILE* b = fopen(init.c_str(), "rb"); if (!b) { fprintf(stderr, "no puedo abrir %s\n", init.c_str()); return 1; }
        int Nf; double Lf; if (fread(&Nf, sizeof Nf, 1, b) != 1 || fread(&Lf, sizeof Lf, 1, b) != 1) return 1;
        if (Nf != N || fabs(Lf - L) > 1e-9) { fprintf(stderr, "%s tiene N=%d L=%g; pide N=%d L=%g\n", init.c_str(), Nf, Lf, N, L); return 1; }
        if (fread(n.data(), sizeof(V), n.size(), b) != n.size()) return 1; fclose(b); }
    double dt = dtf*h*h, tq_fin = thold + 5*tauq, tfin = tq_fin + trampa + relax;
    int pasos = (int)ceil(tfin/dt);
    string nlog = "log_" + tag + ".dat"; FILE* f = fopen(nlog.c_str(), "w");
    fprintf(f, "# nucleacion_persistente.cpp Tp=%g forma=%s Rp=%g Rt=%g k4=%g  N=%d L=%g h=%.4f T0=%g Rc=%g thold=%g tauq=%g theta=%g lam=%g s=%g R0=%g w=%g trampa=%g relax=%g seed=%d init=%s retira=%d\n",
            Tp, forma.c_str(), Rp, Rt, k4, N, L, h, T0, Rc, thold, tauq, theta, lam, s, R0, w, trampa, relax, seed, init.c_str(), retira);
    fprintf(f, "# paso t Tmax lam Q E2 E4 E2+E4 EB fase Q_S R_rms\n");
    printf("tag=%s  N=%d L=%g h=%.4f  dt=%.2e  pasos=%d  (fusión hasta t=%.3f, enfriamiento hasta %.3f, flujo retirado en %.3f, fin %.3f)\n",
           tag.c_str(), N, L, h, dt, pasos, thold, tq_fin, tq_fin + trampa, tfin);
    int nh = omp_get_max_threads(); vector<mt19937_64> rng(nh);
    for (int q = 0; q < nh; q++) rng[q].seed(1000003ULL*seed + 7919ULL*q);
    for (int paso = 0; paso <= pasos; paso++) {
        double t = paso*dt;
        int fase = t < thold ? 1 : (t < tq_fin ? 2 : 3);
        double Tmax = fase == 1 ? T0 : (fase == 2 ? T0*exp(-(t-thold)/tauq) : 0.0);
        double lamE = (fase < 3 || !retira) ? lam : (trampa > 0 ? lam*max(0.0, 1 - (t - tq_fin)/trampa) : 0.0);
        pon_temperatura(Tmax);
        bool registra = paso % cada == 0 || paso == pasos;
        Es e = energia_y_gradiente(n, lamE, true);
        if (registra) { double Q = cargaHopf(n), QS = cargaHopf(n, true), Rr = radio_rms(n);
            fprintf(f, "%d %.5f %.4f %.3f %+.4f %.3f %.3f %.3f %.3f %d %+.4f %.4f\n", paso, t, Tmax, lamE, Q, e.E2, e.E4, e.E2+e.E4, e.EB, fase, QS, Rr); fflush(f);
            if (paso % (cada*20) == 0 || paso == pasos)
                printf("paso %7d t=%.4f fase %d  Tmax=%.3f  λ=%.2f  Q=%+.3f  Q_S=%+.3f  R=%.3f  E2+E4=%.2f\n", paso, t, fase, Tmax, lamE, Q, QS, Rr, e.E2+e.E4); }
        if (paso % cortes == 0 || paso == pasos) guarda_corte(n, paso, t);
        if (paso == pasos) break;
        double amp = sqrt(2*theta*dt/(h*h*h));
        bool ruido = (fase < 3 && Tmax > 1e-3) || Tp > 0;   // el núcleo persistente lleva ruido siempre
        #pragma omp parallel
        { auto& g = rng[omp_get_thread_num()]; normal_distribution<double> nd(0, 1);
          #pragma omp for collapse(2)
          for (int i = 1; i < N-1; i++) for (int j = 1; j < N-1; j++) for (int k = 1; k < N-1; k++) {
              size_t c = I(i,j,k); V gt = G[c] - dot(G[c], n[c])*n[c];
              V mm = n[c] - dt*gt;
              if (ruido) { double T = Tlocal(i, j, k, Tmax);          // ruido proporcional a la temperatura local
                  if (T > 1e-3) mm = mm + (amp*sqrt(min(T/max(T0, Tp), 1.0)))*V{nd(g), nd(g), nd(g)}; }
              n[c] = unit(mm); } }
    }
    fclose(f);
    guarda_bin(n, "n_" + tag + ".bin");
    pon_temperatura(0);                               // queda solo el núcleo persistente (si Tp > 0)
    double Q = cargaHopf(n), QS = cargaHopf(n, true), Rr = radio_rms(n); Es e = energia_y_gradiente(n, 0, false);
    printf("FINAL %s: Q = %+.3f  Q_S = %+.3f  R_rms = %.3f  E2+E4 = %.2f  E2/E4 = %.3f\n",
           tag.c_str(), Q, QS, Rr, e.E2+e.E4, e.E4 > 0 ? e.E2/e.E4 : 0.0);
    FILE* r = fopen("resumen.dat", "a");
    fprintf(r, "%-22s lam=%g s=%+g retira=%d T0=%g Tp=%g forma=%s Rp=%g k4=%g seed=%d N=%d L=%g init=%s  Q=%+.3f  Q_S=%+.3f  R_rms=%.3f  E2+E4=%.2f  E2/E4=%.3f\n",
            tag.c_str(), lam, s, retira, T0, Tp, forma.c_str(), Rp, k4, seed, N, L, init.c_str(), Q, QS, Rr, e.E2+e.E4, e.E4 > 0 ? e.E2/e.E4 : 0.0);
    fclose(r);
}
