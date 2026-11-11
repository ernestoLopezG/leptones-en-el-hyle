# Leptones en el hyle — programas

Código que acompaña al preprint «Leptones en el hyle» (Ernesto López). Todo el código lo ha escrito Claude (Anthropic) y lo ha revisado y ejecutado el autor. Las cabeceras de cada programa explican la ecuación que resuelve, el método y el formato de salida.

DOI: [10.5281/zenodo.23545994](https://doi.org/10.5281/zenodo.23545994) (todas las versiones) · Licencia: MIT (fichero `LICENSE`).

## Requisitos (Debian / Ubuntu)

```
sudo apt install g++ python3-numpy python3-scipy python3-matplotlib python3-skimage libfftw3-dev
```

Compilar cualquiera de los programas en C++ (todos son un único fichero y usan OpenMP):

```
g++ -O3 -march=native -fopenmp programa.cpp -o programa
```

Para elegir los núcleos: `OMP_NUM_THREADS=4 ./programa ...`

## Unidades

Todo es adimensional.

- **Carpetas 1 a 4:** ħ = c = m = 1.
- **Carpetas 5, 7 y 8:** las longitudes están en unidades de ℓ = √(κ₄/κ₂), k en unidades de 1/ℓ y las energías en unidades de √(κ₂κ₄). κ₂ es la rigidez del director y κ₄ la del término de cuarto orden. El valor de ℓ en metros no está fijado; es el punto 5 de «Lo abierto» del preprint (sección 6).
- **Carpeta 6:** la ecuación de Mathieu es adimensional.

## Qué hay en cada carpeta y a qué parte del preprint corresponde

| Carpeta | Programas | Preprint |
|---|---|---|
| `1_reposo` | `reposo2.cpp`: leptón en reposo con simetría esférica (programa de referencia de la fase escalar) | Sección 5, contexto |
| `2_umbral` | `umbral.cpp` (estática con umbral inferior tipo Bingham), `evolucion.cpp` (frecuencia libre), `anclada.cpp` (frecuencia anclada), `derrick.cpp` (sin umbral), `cargafija.cpp` y `arrastre.cpp` (pruebas auxiliares) | Sección 5, filas 1 y 2 |
| `3_giro` | `giro.cpp`, `giro_anc.cpp`: envolvente que gira en un plano de las extendidas | Sección 5, fila 4 |
| `4_cuartico` | `cuartico.cpp`: solitón de umbral con término de cuarto orden, sin topología | Sección 5, fila 3 |
| `5_hopfion` | `hopfion.cpp`, `hopfion_quiral.cpp`, `hopfion_P.cpp`, los análisis en Python, `perturba.py` (punto de silla) y `convergencia.sh` | Sección 4; sección 5, filas 5 y 6; errores corregidos |
| `7_nucleacion` | `nucleacion.cpp` (formación del nudo por fusión y flujo helicoidal), `nucleacion_persistente.cpp` y `nucleacion_pulsos.cpp` (el nudo frente a la fusión), `torsion_flujo.py` y los barridos | Sección 4, «Formación del nudo» y «El nudo frente a la fusión»; sección 5 |
| `8_tamano` | `hopfion_compacto.cpp` (dirección periódica que se encoge) y `hopfion_capas.cpp` (lámina de esméctico con capas; necesita FFTW) | Sección 5, vías descartadas sobre el tamaño del nudo |
| `6_mathieu` | `mathieu_q.py`: funciones angulares de Mathieu con q negativo (método de Hill) | Sección 2, figura de los pozos |
| `preprint` | `leptones_en_el_hyle.tex` (fuente LaTeX), `leptones_en_el_hyle.pdf` y las figuras en `figuras/` | El texto completo |

## Compilar el preprint

El texto está en `preprint/leptones_en_el_hyle.tex` y se compila con XeLaTeX (dos pasadas, para las referencias cruzadas):

```
sudo apt install texlive-xetex texlive-latex-extra texlive-lang-spanish texlive-science
cd preprint
xelatex leptones_en_el_hyle.tex
xelatex leptones_en_el_hyle.tex
```

## Versión 1.1 (noviembre de 2026)

La versión 1.1 retira del preprint el material del hopfión con giro k ≠ 0 en η (masa de Kaluza–Klein, cola de Yukawa, forma de la cola): suponía un giro estático, incompatible con el postulado 2. Los programas `hopfion_eta.cpp`, `cola_angular.py` y `figura_cola_dipolar.py` y las órdenes con k ≠ 0 de esta guía se mantienen como parte de la versión 1.0. Añade la tabla de convergencia, la prueba de punto de silla, la formación del nudo por nucleación (carpeta 7) y las vías descartadas sobre su tamaño (carpeta 8). En el preprint 1.1 la numeración de secciones cambia: 3 es la nueva sección de la fuerza eléctrica, 4 la estructura en las extendidas, 5 lo descartado, 6 lo abierto y 7 las conclusiones.

## Convergencia y punto de silla (sección 4)

En la carpeta `5_hopfion`:

```
bash convergencia.sh
```

| Malla | Caja | h | E/(32√2π²) | E2/E4 | Q |
|---|---|---|---|---|---|
| 61³ | 5,4 | 0,090 | 1,280 | 0,66 | −0,961 |
| 81³ | 7,2 | 0,090 | 1,241 | 0,78 | −0,968 |
| 101³ | 9,0 | 0,090 | 1,224 | 0,86 | −0,971 |
| 121³ | 9,0 | 0,075 | 1,232 | 0,87 | −0,980 |

Resultado: E/(32√2π²) = 1,23 ± 0,02 (literatura ≈ 1,22). Punto de silla: el nudo deformado (+46 % de energía) vuelve a E/(32√2π²) = 1,27995 con simetría axial.

## Formación del nudo y fusión (sección 4; carpeta `7_nucleacion`)

```
bash barrido.sh rapido          # formación desde el vacío (barrido.sh … barrido5.sh: barridos de parámetros)
bash barrido_persistente.sh     # núcleo fundido permanente, forzado ½ frente a 3/2 (unas 2,5 h con 4 núcleos)
bash barrido_pulsos.sh          # fusión pulsada (unos 30 min)
```

Resultados esperados: el signo de la helicidad del flujo elige el signo de Q y la carga crece con la torsión; una fusión permanente sobre el nudo lo desata (incluso con κ₄ y sin ruido); con pulsos breves sobrevive; en medio ordenado el nudo Q = 1 no cambia de carga bajo forzado alineado. Las cabeceras de cada programa explican los parámetros.

## Tamaño del nudo (sección 5; carpeta `8_tamano`)

```
bash barrido_compacto.sh        # al encoger la dirección periódica el nudo se aplasta pero no se reduce en las extendidas
bash barrido_capas.sh           # el radio no sigue a √(λΔ): R/R_control ≈ 0,90
```

## Reproducir los resultados de la sección 4 (versión 1.0: sección 3)

Todo se hace en la carpeta `5_hopfion`. Los tiempos son orientativos para un i5 de 4 núcleos.

**Hopfión de referencia** (malla 150³, caja de lado 12; muchas horas, se puede cortar y continuar):

```
OMP_NUM_THREADS=4 ./hopfion 150 12 1.0 30000
python3 hopfq.py n_final_N150.bin
```

Resultado esperado: E = 544,35, E/(32√2π²) = 1,219 (literatura: ≈ 1,22), Q = −0,979, E2/E4 = 0,965.

**Giro en η** (versión 1.0; retirado del preprint 1.1):

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

Da el reparto de la energía del flujo: 66 % poloidal, 34 % toroidal (velocidad en el gauge de Coulomb).

**Figura de la sección 4:**

```
python3 figura_hopfion.py n_final_N150.bin
```

(La figura con giro en η, `n_eta_N150_K1.bin`, es de la versión 1.0.)

**Tamaño y cola** (radio del anillo y ajuste de Yukawa k_ef; la cola es material de la versión 1.0):

```
python3 medir_tamano.py n_final_N150.bin n_eta_N150_K0.5.bin n_eta_N150_K1.bin
```

| k | R_anillo | k_ef |
|---|---|---|
| 0 | 1,184 | — |
| 0,5 | 1,094 | — |
| 1 | 0,835 | 1,064 |

Para k = 0 y k = 0,5 el valor de k_ef lo contamina la caja (ver «Controles»).

**Forma angular de la cola** (versión 1.0; perfil dipolar sin²θ alrededor del eje del nudo):

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

**Precesión del eje bajo un campo externo** (figura de la sección 6, punto 3):

```
python3 precesion_hopfion.py n_final_N150.bin 30    # el último número es el ángulo del cono en grados
```

**Factor g orbital del director** (el cálculo corregido, que da g = 1):

```
python3 gfactor.py n_final_N150.bin
```

## Funciones de Mathieu (sección 2)

En la carpeta `6_mathieu`:

```
python3 mathieu_q.py    # → mathieu_q_negativo.png, con q = −10, −100 y −1000
```

## Reproducir la sección 5 (vías descartadas)

Los argumentos de estas órdenes son ejemplos; las cabeceras explican cada parámetro.

```
./umbral 0.1 0.3 1                      # soluciones estáticas con umbral uc
./evolucion 1 0.01 200                  # frecuencia libre: cae a un estado con autoenlace (ω ≈ 0,88)
./anclada 1 0.01 200                    # frecuencia anclada
./derrick 2 0.01 200                    # sin umbral (teorema de Derrick)
./cuartico 1 1 0.3 0 200                # término de cuarto orden sin topología
./giro 1 12 2 200 1                     # giro en un plano (ANCLADA=1 ./giro_anc ... para la versión anclada)
K=1 C4=0 Q0=0.8 ./hopfion_quiral 100 8 1.0 20000    # quiralidad en lugar del cuarto orden: el nudo se deshace
K=1 S=0 C4=1 ./hopfion_P 81 6 0.9 2000              # control: reproduce el nudo con giro (Q ≈ −0,97)
K=1 S=90 C4=0 ./hopfion_P 81 6 0.9 8000 n_P_N81_K1_S0_C1.bin   # giro en el tiempo sin cuarto orden: ω sube a ≈ 1,7 > k y el nudo se deshace (Q → 0)
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
- **Paso de tiempo:** el flujo de gradiente es explícito. Con nudos pequeños (k grande), el término de cuarto orden exige pasos muy cortos. `hopfion_eta`, `hopfion_quiral` y `hopfion_P` reducen el paso solos cuando la energía sube y lo avisan en pantalla. `DTAU=...` fija el paso inicial en unidades de h².
- **k ≥ 2 (versión 1.0):** el nudo se contrae por debajo de lo que resuelven mallas con h ≈ 0,06 y se pierde. Era una limitación declarada en la versión 1.0 del preprint; el material con k ≠ 0 se ha retirado.
- **Cola en la caja de lado 12:** con k = 0 y k = 0,5 la pared contamina el ajuste de la cola. Por ejemplo, k = 0 da un k_ef = 0,37 espurio.
- **Ficheros `.dat` antiguos:** la columna de energía normalizada es ahora E/(32√2π²), en la normalización de Sutcliffe (2007). Los `.dat` generados antes de esta versión traían E/504.

## Referencias

- P. Sutcliffe, «Knots in the Skyrme–Faddeev model», Proc. R. Soc. A 463, 3001 (2007), arXiv:0705.1468.
- R. A. Battye y P. M. Sutcliffe, PRL 81, 4798 (1998).

## 9_burbuja — el nudo en una región ordenada finita (v1.1, revisión tras el arbitraje)

`hopfion_burbuja.cpp` relaja un hopfión dentro de una región esférica ordenada de radio Rb (borde de anchura w, orden S(r) que cae a cero), sin borde fijo y con la rigidez pesada por S² y S⁴; fuera, fondo uniforme o de direcciones aleatorias. `hopfq_burbuja.py` calcula Q con ventana. `barrido_burbuja.sh` hace el barrido del preprint (unas 5 h con 4 núcleos).

Resultado (h ≈ 0,09, w = 0,4): con Rb = 2 y 3 el nudo se escapa; con Rb = 4 está escapándose; con Rb = 6 y 8 sobrevive (E/(32√2π²) ≈ 1,16–1,19, Q ≈ −0,97). El fondo aleatorio da lo mismo que el uniforme.

## 10_polaron — la circulación repartida entre la onda y el patrón (v1.1, segundo arbitraje)

Juguete unidimensional en el anillo η ∈ [0, 2π). Una onda antiperiódica (Schrödinger, ħ = 1, radio 1) siente un potencial −λφ, y el medio φ, con inercia M, responde a la intensidad: M φ_tt = −κφ + Dφ_ηη + λ|Ψ|². Las dos ecuaciones salen de una sola energía, así que E y el momento total P = P_onda + P_medio, con P_medio = −M∫φ_t φ_η, se conservan. Se trabaja en el sector de dos zonas (componentes k = ½ + 2n).

- `polaron_anillo.py` calcula el estado estático (patrón quieto) y el que deriva a velocidad angular v con P = ½. Usa la iteración autoconsistente en el sistema que deriva y una bisección en v. κ = 1 y D = 0,2. Escribe `resultado.txt` (unos 5 min).
- `evolucion.py lam M` evoluciona en el tiempo el estado que deriva (paso partido para la onda, Verlet para el medio). Comprueba que P se conserva y que el patrón gira a la v prevista (`evolucion_M100.txt`, unos 20 s).

Resultado: el patrón estático deja a la onda solo parte del momento (0,355 con λ² = 10; 0,232 con λ² = 12; 0,104 con λ² = 16). El estado con P = ½ existe con el patrón intacto. El medio lleva el resto, la deriva es v ∝ 1/M y el coste de energía de la deriva también es ∝ 1/M. Con M = 100 y λ² = 10: v = 0,0025 y ΔE = 1,8×10⁻⁴. En la evolución, P se conserva y el patrón gira a la v prevista hasta la quinta cifra.

Límite del juguete: en una teoría clásica P puede tomar cualquier valor. El ½ se impone aquí como objetivo; su cuantización exacta sale del argumento de simetría del preprint (antiperiodicidad del conjunto onda más patrón), que es cuántico.

`10_polaron/solapamientos.py` calcula el solapamiento |⟨Ψ(η)|Ψ(η − π/2)⟩| de dos ondas de Mathieu de orden ½ con las líneas de zonas de mayor índice perpendiculares (imagen de dos electrones de espín opuesto). Da 0,83 con q = −0,64, 0,69 con q = −1 y 0,05 con q = −10. Con el |q| moderado que se espera para el electrón, la separación geométrica es parcial.

`preprint/figuras_src/fig00_ensamblaje.py` dibuja la figura de ensamblaje del electrón (espacio ordinario y disco compacto en tres puntos).
