# Leptones en el hyle — programas

Código que acompaña al preprint «Leptones en el hyle» (Ernesto López). Los programas los escribió Claude (Anthropic) en colaboración con el autor; en las cabeceras van marcados con †. Las cabeceras de cada programa explican la ecuación que resuelve, el método y el formato de salida.

DOI: [10.5281/zenodo.23545995](https://doi.org/10.5281/zenodo.23545995) · Licencia: MIT (fichero `LICENSE`).

## Requisitos (Debian / Ubuntu)

```
sudo apt install g++ python3-numpy python3-scipy python3-matplotlib python3-skimage
```

Compilar cualquiera de los programas en C++ (todos son un único fichero y usan OpenMP):

```
g++ -O3 -march=native -fopenmp programa.cpp -o programa
```

Para elegir los núcleos: `OMP_NUM_THREADS=4 ./programa ...`

## Unidades

Todo es adimensional.

- **Carpetas 1 a 4:** ħ = c = m = 1.
- **Carpeta 5:** las longitudes están en unidades de ℓ = √(κ₄/κ₂), k en unidades de 1/ℓ y las energías en unidades de √(κ₂κ₄). κ₂ es la rigidez del director y κ₄ la del término de cuarto orden. El valor de ℓ en metros no está fijado; es el punto 5 de «Lo abierto» del preprint.

## Qué hay en cada carpeta y a qué parte del preprint corresponde

| Carpeta | Programas | Preprint |
|---|---|---|
| `1_reposo` | `reposo2.cpp`: leptón en reposo con simetría esférica (programa de referencia de la fase escalar) | Sección 4, contexto |
| `2_umbral` | `umbral.cpp` (estática con umbral inferior tipo Bingham), `evolucion.cpp` (frecuencia libre), `anclada.cpp` (frecuencia anclada), `derrick.cpp` (sin umbral), `cargafija.cpp` y `arrastre.cpp` (pruebas auxiliares) | Sección 4, filas 1 y 2 |
| `3_giro` | `giro.cpp`, `giro_anc.cpp`: envolvente que gira en un plano de las extendidas | Sección 4, fila 4 |
| `4_cuartico` | `cuartico.cpp`: solitón de umbral con término de cuarto orden, sin topología | Sección 4, fila 3 |
| `5_hopfion` | `hopfion.cpp`, `hopfion_eta.cpp`, `hopfion_quiral.cpp` y los análisis en Python | Sección 3; sección 4, fila 5; errores corregidos |

## Reproducir los resultados de la sección 3

Todo se hace en la carpeta `5_hopfion`. Los tiempos son orientativos para un i5 de 4 núcleos.

**Hopfión de referencia** (malla 150³, caja de lado 12; muchas horas, se puede cortar y continuar):

```
OMP_NUM_THREADS=4 ./hopfion 150 12 1.0 30000
python3 hopfq.py n_final_N150.bin
```

Resultado esperado: E = 544,35, E/(32√2π²) = 1,219 (literatura: ≈ 1,21), Q = −0,979, E2/E4 = 0,965.

**Giro en η**:

```
K=0.5 OMP_NUM_THREADS=4 ./hopfion_eta 150 12 1.0 30000
K=1   OMP_NUM_THREADS=4 ./hopfion_eta 150 12 1.0 30000
```

| k | E | Q |
|---|---|---|
| 0,5 | 586,68 | −0,973 |
| 1 | 675,38 | −0,959 |

**Flujo de Mermin–Ho:**

```
python3 flujo.py n_final_N150.bin
```

Da el reparto de la energía del flujo: 66 % poloidal, 34 % toroidal.

**Figuras de la sección 3:**

```
python3 figura_hopfion.py n_final_N150.bin
python3 figura_hopfion.py n_eta_N150_K1.bin
```

**Tamaño y cola** (radio del anillo y ajuste de Yukawa k_ef):

```
python3 medir_tamano.py n_final_N150.bin n_eta_N150_K0.5.bin n_eta_N150_K1.bin
```

| k | R_anillo | k_ef |
|---|---|---|
| 0 | 1,184 | — |
| 0,5 | 1,094 | — |
| 1 | 0,835 | 1,064 |

Para k = 0 y k = 0,5 el valor de k_ef lo contamina la caja (ver «Controles»).

**Forma angular de la cola** (perfil dipolar sin²θ alrededor del eje del nudo):

```
python3 cola_angular.py n_final_N150.bin n_eta_N150_K0.1.bin
RIN=2 python3 cola_angular.py n_eta_N150_K1.bin
```

Resultado esperado: a2/a0 ≈ −1,00 desde r ≈ 3 con k = 0 y k = 0,1; con k = 1, a2/a0 ≈ −1,08 y a4/a0 ≈ +0,09, un octupolo de un 0,8 % que decae como k₃(kr)/k₁(kr).

**Reparto de la energía del flujo** (poloidal frente a toroidal, en el mapa de Hopf analítico y en los .bin):

```
python3 reparto_flujo.py n_final_N150.bin n_eta_N150_K0.5.bin n_eta_N150_K1.bin
```

Resultado: el mapa analítico tiende a 0,50 / 0,50 al afinar la malla; el hopfión da 0,658 de poloidal con k = 0, 0,667 con k = 0,5 y 0,679 con k = 1. El reparto depende de la forma del mínimo de energía, no de la topología.

**Figuras ilustrativas** (campos analíticos, no la simulación):

```
python3 figura_cola_dipolar.py    # cola de Yukawa dipolar alrededor de un nudo pequeño
python3 figura_espin_flujo.py     # líneas de flujo para espín +½, −½ e imagen especular, con su enlace
```

**Factor g orbital del director** (el cálculo corregido, que da g = 1):

```
python3 gfactor.py n_final_N150.bin
```

## Reproducir la sección 4 (vías descartadas)

Los argumentos de estas órdenes son ejemplos; las cabeceras explican cada parámetro.

```
./umbral 0.1 0.3 1                      # soluciones estáticas con umbral uc
./evolucion 1 0.01 200                  # frecuencia libre: cae a un estado con autoenlace (ω ≈ 0,88)
./anclada 1 0.01 200                    # frecuencia anclada
./derrick 2 0.01 200                    # sin umbral (teorema de Derrick)
./cuartico 1 1 0.3 0 200                # término de cuarto orden sin topología
./giro 1 12 2 200 1                     # giro en un plano (ANCLADA=1 ./giro_anc ... para la versión anclada)
K=1 C4=0 Q0=0.8 ./hopfion_quiral 100 8 1.0 20000    # quiralidad en lugar del cuarto orden: el nudo se deshace
```

## Formato de los ficheros `.bin`

Son ficheros binarios con este contenido, en este orden:
- un `int32` con N;
- un `double` con L;
- N³ vectores (n_x, n_y, n_z) en `double`, con índice (i, j, k) y k el más rápido.

La malla es x = −L/2 + h·i, con h = L/(N − 1). El borde está fijo a n = (0, 0, 1).

Los tres programas de hopfión guardan el `.bin` cada 2000 iteraciones. Para continuar un cálculo cortado, se añade el `.bin` como último argumento, con los mismos N y L.

## Controles y limitaciones conocidas

- **Virial:** en el mínimo sin giro, E2/E4 → 1. Con giro, E2 − E4 + 3·Et2 + Et4 → 0, que es la última columna de `hopfion_eta`.
- **Cota de energía:** los términos del giro son positivos. Por eso, mientras Q = 1, la energía con giro nunca puede bajar de la de k = 0 (544). Si baja, el nudo se ha deshecho en la malla.
- **Malla:** con h = L/(N − 1) mayor que unos 0,09 el nudo se deshace, incluso con k = 0.
- **Paso de tiempo:** el flujo de gradiente es explícito. Con nudos pequeños (k grande), el término de cuarto orden exige pasos muy cortos. `hopfion_eta` y `hopfion_quiral` reducen el paso solos cuando la energía sube y lo avisan en pantalla. `DTAU=...` fija el paso inicial en unidades de h².
- **k ≥ 2:** el nudo se contrae por debajo de lo que resuelven mallas con h ≈ 0,06 y se pierde. Es una limitación declarada en el preprint.
- **Cola en la caja de lado 12:** con k = 0 y k = 0,5 la pared contamina el ajuste de la cola. Por ejemplo, k = 0 da un k_ef = 0,37 espurio.
- **Ficheros `.dat` antiguos:** la columna de energía normalizada es ahora E/(32√2π²), en la normalización de Sutcliffe (2007). Los `.dat` generados antes de esta versión traían E/504.

## Referencias

- P. Sutcliffe, «Knots in the Skyrme–Faddeev model», Proc. R. Soc. A 463, 3001 (2007), arXiv:0705.1468.
- R. A. Battye y P. M. Sutcliffe, PRL 81, 4798 (1998).
