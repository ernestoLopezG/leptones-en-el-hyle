# Informe de arbitraje — «Leptones en el hyle» (versión 1.1)

## (a) Resumen

El trabajo propone una imagen mecánica del electrón: una pulsación sin masa de un medio continuo de tipo esméctico en 5+1 dimensiones. La pulsación circula a velocidad *c* por una dirección compacta angular η, en la interfaz entre una fase A y una fase C\*. La masa es la frecuencia de corte de esa guía de ondas. La carga es la circulación en la interfaz. *g* = 2 sale de postular que el marco local está anclado a las capas (conexión de Weitzenböck). La ley de Coulomb se obtiene de un arrastre del medio a lo largo de η. En las tres dimensiones extendidas, una textura del director con carga de Hopf 1 (un hopfión de Faddeev–Skyrme) marca el centro y el eje del espín. Esa parte numérica es la más sólida: es reproducible y concuerda con la literatura.

El texto es honesto sobre lo que no deduce (α, *g*−2, las masas) y etiqueta casi todo. Aun así, tiene un problema conceptual de fondo: identifica con el espín una magnitud que, en su propia reducción de Kaluza–Klein, es la carga eléctrica, es decir, el número cuántico de la circulación en η. De esa identificación salen el «ħ/2 cinemático», la exclusión de los «estados de espín 3/2» y la regla k_director = 2k_espinor, y ninguna de las tres se sostiene tal como está escrita. Además, he comprobado numéricamente que las nucleaciones (q ≠ 0), en las que descansan la lectura del ángulo del espín y la imagen de Pauli, reducen de forma exponencial la circulación media en η, y por tanto la carga y el momento angular que el modelo atribuye a la pulsación. Hay también varias etiquetas [D] y [N] que no cumplen las definiciones del propio trabajo, y una ecuación de ondas en el plano compacto que nunca se escribe y con la que no se pueden comprobar las afirmaciones sobre perfiles y masas.

---

## (b) Problemas graves

### G1. Se confunde la carga (momento en η) con el espín (rotaciones en 3D)
**Dónde:** postulado 2, «Orden ½», Pasos 2–5 de §2.1, punto 5 de la lista de §2.1, punto 9(iii) y punto 10 de §6, Conclusiones.

En Kaluza–Klein, el número cuántico de la circulación en la dirección compacta, ν (o k = ν/ξ₀), es la **carga**: es el generador de traslaciones en η, una simetría U(1). El propio texto lo usa así en el Paso 4 («acoplamiento mínimo de una carga q ∝ k»). Ese número no se transforma bajo las rotaciones del espacio de 3 dimensiones y no puede ser, a la vez, el momento angular de espín. Consecuencias concretas:

1. **El «ħ/2» del postulado 2** (m ξ₀ c = ħ/2) es ħν, el momento en η de un modo con ν = ½: es decir, la carga en otras unidades. En la reducción del Paso 5, el espín ½ sale de que Ψ es un espinor de Dirac («Y como Ψ es un espinor, el espín es ½»), no de ν. El texto lo reconoce a medias («Este momento angular está en el plano compacto; su dirección en las extendidas la marca el eje del nudo [H]»), pero conclusiones y resumen presentan la identidad como una explicación del espín («sin ninguna superficie que supere a c, lo que permite superar las contradicciones originales de 1925»). Lo que de verdad resuelve el problema de 1925 en el texto es el temblor de Dirac (Huang, Hestenes), que ya es estándar.
2. **«Orden 3/2 = espín 3/2» (§6, punto 10; Conclusiones, «Tumbarían el núcleo…»).** Con la misma reducción del Paso 3, un modo ν = 3/2 da una ecuación de Dirac con masa 3m y carga ∝ 3/2: es un **fermión de espín ½** con carga 3e, no un estado de espín 3/2. La torre de Kaluza–Klein de un espinor en 5D está formada solo por fermiones de Dirac de espín ½. El problema real que hay que resolver es la ausencia de copias de carga 3e y masa ≈ 1,5 MeV, que serían muy fáciles de producir. El texto lo dice en un sitio («tendría carga 3e») y en otro habla de «un estado del electrón de espín 3/2 o excitado, con su misma carga». Es contradictorio.
3. **«Ψ es un espinor de cuatro componentes: dos sentidos de circulación (carga) por dos orientaciones (espín)»** (punto 5 de la lista). Es incorrecto: el Paso 2 fija el signo de k, y χ sigue teniendo cuatro componentes. Son las cuatro de un espinor de Dirac en 3+1 dimensiones (dos de espín y dos de energía positiva y negativa). El modo −k es el conjugado de carga, no «las otras dos componentes».
4. **El argumento del orden ½** (§2, «Orden ½») hace girar el marco en el plano (ξ, η), pero el Paso 1 usa un espinor de 4+1 dimensiones, sobre el que no actúan las rotaciones en ese plano. El texto reconoce que haría falta un espinor de 5+1 dimensiones (8 componentes), y no lo usa. Tal como está, el orden ½ es una condición de contorno elegida (condiciones de Neveu–Schwarz), lo cual es legítimo, pero no está «motivado por los postulados 1, 3, 4 y 6».

**Propuesta:** separar sin ambigüedad tres cosas: (i) ν = ½ como condición antiperiódica (estructura de espín), etiquetada [H]; (ii) la carga ∝ ⟨−i∂_η⟩; y (iii) el espín ½ como propiedad del espinor de Dirac en 3+1 dimensiones. Retirar «ħ/2 = m ξ₀ c» como explicación del espín, o reformularlo como «el temblor de Dirac tiene radio ξ₀ y velocidad c», que es estándar. Sustituir «espín 3/2» por «modo ν = 3/2: fermión de carga 3e y masa 3m» y discutir con datos por qué no se observa.

### G2. Las nucleaciones suprimen la circulación, y con ella la carga y el «ħ/2»
**Dónde:** §2, «Origen de la ecuación de Mathieu», «Por qué el orden ½ puede circular», figura 3, «Dos electrones por orbital», Paso 4.

Para un estado de Floquet Ψ = e^{iνη}P(η) de la ecuación de Mathieu, el momento medio en η cumple, por Hellmann–Feynman, ⟨−i∂_η⟩ = ½ ∂a/∂ν. En el centro de la banda (ν = ½) es proporcional a la pendiente de la banda, que con |q| grande es exponencialmente pequeña. Lo he calculado por el método de Hill (201 armónicos) y coincide con la derivada numérica de a(ν):

| q | a_{1/2} | ⟨−i∂_η⟩ (valor sin potencial: 0,5) |
|---|---|---|
| 0 | 0,25 | 0,500 |
| −0,6393 | 1×10⁻⁵ | 0,366 |
| −1 | −0,307 | 0,262 |
| −3 | −2,81 | 3,8×10⁻² |
| −10 | −13,94 | 3,4×10⁻⁴ |
| −30 | −49,3 | 7,5×10⁻⁸ |

(De paso: confirma que q = −0,6393 da a ≈ 0, como afirma el texto para el tercer modo, y que la semianchura de la primera banda con q = −10 es 2,1×10⁻⁴.)

Cuando se sustituye ∂_η por ik en el acoplamiento ∂_μ − κA_μ∂_η, la carga efectiva en 3+1 dimensiones es κ⟨−i∂_η⟩, no κν. Del mismo modo, el momento angular en el plano compacto es ħ⟨−i∂_η⟩, no ħ/2. Con los valores de |q| que usa el texto para la imagen del espín y de Pauli (q = −10 a −1000), la «circulación de fase» es casi nula. La afirmación «el caudal j es constante en toda la vuelta» es cierta, pero ese caudal constante es proporcional a ⟨−i∂_η⟩, que es diminuto. La universalidad de la carga exige además que todas las partículas cargadas tengan exactamente el mismo q, lo que no se argumenta. Por último, a_{1/2} < 0 para |q| > 0,64: si a hace el papel de (k_c ξ₀)², el corte se vuelve imaginario. Hay que explicar qué relación hay entre a y la masa.

**Propuesta:** calcular la carga y el momento angular como valores esperados con el modo completo, que incluye la contribución del propio patrón de nucleaciones, porque al ser autogenerado se mueve con la onda y puede llevar su propio momento. Mostrar que la suma total da exactamente e y ħ/2 con independencia de q. Si no sale, el electrón tiene que tener |q| ≲ 1, y entonces la imagen de «dos nucleaciones casi aisladas» y los solapamientos de 10⁻⁵–10⁻¹⁵ del punto 7 de §2.1 dejan de valer.

### G3. La ecuación en el plano compacto no está escrita; perfiles y masas no se pueden comprobar
**Dónde:** «Separación de variables», tabla de tipos de §2.0, «Masas de los neutrinos», punto 2 de la lista de §2.1, figura 1.

- Si el operador del plano compacto es el laplaciano bidimensional ordinario, un perfil I_ν(κξ)e^{iνη} cumple ∇²Φ = +κ²Φ **exactamente**: el término ν²/ξ² ya está dentro de la ecuación de Bessel. El corte sería k_c² = −κ² < 0, es decir, taquiónico. Para que las Bessel modificadas den masa real, el medio tiene que tener un operador anisótropo o un potencial radial, y eso es justo lo que falta. La descomposición «ω² = c²k_arco² + c²k_ξ²» del punto 2 es una lectura local (tipo WKB), no la solución del problema de valores propios.
- El texto dice tres cosas incompatibles sobre ξ₀: «la relación ξ₀ = ħ/(2mc) vale tal cual para el electrón […] y fija el radio de la interfaz: ξ_interfaz ≈ 1,93×10⁻¹³ m» (Paso 3, glosario, figura 1); «el radio físico de la interfaz, ξ_i, es algo menor que ξ₀» (punto 2); y «la relación ξ₀ = ħ/(2mc) es el límite q → 0» («La masa con q ≠ 0»).
- Las masas de los neutrinos, que son el único resultado cuantitativo del sector compacto y la base de la calibración de R_u, no se dan como números («del orden de las centésimas de eV»), no se muestra el cálculo y llevan la etiqueta [N] con la nota «cálculo no incluido en el repositorio». La condición de cuantización («el baricentro en ξ/ξ₀ = 1») es *ad hoc* y no tiene justificación física.

**Propuesta:** escribir la ecuación de ondas del hyle en el plano compacto (con la anisotropía, el potencial o las condiciones de contorno de interfaz) y resolver al menos el caso de los neutrinos tipo I, dando m₁, m₂ y R_u = 5,1 µm con el código en el repositorio. Unificar la definición de ξ₀ (definida por la masa) y la de ξ_i (radio físico) en todo el texto.

### G4. La topología del plano compacto: disco para el espín, corona para el tercer neutrino
**Dónde:** «Orden ½» de §2 frente a §2.0 («posible gracias al hueco central»), Conclusiones.

El argumento del orden ½ necesita que cada capa «rodee un disco» («Los modos enteros solo serían posibles si la circunferencia no rodeara nada, como en un toro sin centro»). Pero el modo IIbis solo existe porque «no atraviesa el centro: el hueco central, el núcleo fluido». Y en las conclusiones el núcleo fluido «es también lo que permite los hadrones». Si el núcleo fluido actúa como un agujero para las ondas, el dominio es topológicamente una corona circular, que admite tanto la estructura de espín periódica como la antiperiódica, y el argumento del disco se cae. Si no actúa como agujero, la solución IIbis «no sería válida en un disco lleno», según el propio texto. Además, el postulado 6 dice que transportar el marco por la capa «no lo gira», y el argumento del orden ½ dice que «al dar una vuelta gira 360° respecto a un marco fijo». Las dos frases pueden ser ciertas a la vez, porque se refieren a marcos distintos, pero el texto debe decir respecto a qué marco se define la (anti)periodicidad del espinor.

**Propuesta:** decidir si la estructura de espín se extiende por el núcleo fluido, y declarar el orden ½ como postulado independiente [H] sin presentarlo como consecuencia geométrica.

### G5. La ley de Coulomb: el argumento de «sin término de masa» confunde un desplazamiento con una velocidad
**Dónde:** §3, «La ecuación del campo: por qué 1/r»; tabla final de §3; resumen.

El texto argumenta: «desplazar su origen de manera uniforme no cambia nada físico… Por tanto un arrastre uniforme no cuesta energía». La simetría prohíbe términos en el desplazamiento λ, pero u es una **velocidad** del medio (∂_tλ). En un medio con inercia, que es el caso porque se postula superfluido como la componente de Landau, un flujo uniforme tiene energía cinética ½ρu². Ese término es precisamente un término de masa, y con él la ecuación sería ∇²u = (ρ/K_u)u y el potencial de Yukawa, no 1/r. En Kaluza–Klein, A₀ → A₀ − ∂_tλ es una transformación de coordenadas sin coste. Aquí, en cambio, el hyle tiene un sistema en reposo físico («el tiempo es […] el del sistema en reposo del hyle»), y una rotación uniforme del medio en η es, en principio, física. Hace falta un argumento explícito de por qué su energía cinética no cuenta. El propio texto usa este mecanismo para la neutrínica («un arrastre uniforme sí cuesta energía»), de modo que tiene que explicar por qué en la interfaz no ocurre.

Otros problemas de la misma sección:
- La simetría de la fuerza (tercera ley) exige que la amplitud de la fuente A sea proporcional a la carga de la fuente. No se deduce: solo se iguala para dos electrones.
- El signo se etiqueta «D por analogía hidrodinámica», cuando el propio texto admite que «que el arrastre en las extendidas […] conserve el signo […] no está calculado [A]». Con una parte [A], el conjunto no puede ser [D]. Debería ser [H].
- El resumen dice que «se deducen la ley de Coulomb y su signo». Con α igualada, 1/r dependiente del punto anterior y el signo por analogía, la palabra «deducen» es excesiva.

**Propuesta:** escribir la energía del medio (cinética más gradiente) para el campo de arrastre y mostrar qué la deja sin masa en la interfaz. Si no se consigue, etiquetar 1/r como [H].

### G6. El director como bilineal del espinor no puede girar con la fase de carga
**Dónde:** postulado 7 («como es un bilineal del espinor, su número de onda es el doble»), §6 punto 3, punto 5(iii).

El mapa de Hopf **n** = z†**σ**z es invariante ante z → e^{iα}z. Precisamente por eso es una fibración: elimina la fase global. Si avanzar en η es la fase global del espinor, que es lo que el Paso 4 identifica con la carga y con el gauge, el director no gira en absoluto, y la relación k_director = 2k_espinor es falsa. Solo giraría si el avance en η fuera una fase **relativa** entre las dos componentes, es decir, una rotación de espín alrededor del eje. Eso devuelve el texto a la confusión de G1: o la fase de η es la carga, o es una rotación de espín, pero no ambas. De ahí dependen también «el giro con la pulsación», la cota (iii) del punto 5 y la cancelación de los términos de giro que justifica usar el hopfión puro en 3D.

**Propuesta:** escribir explícitamente el espinor y el bilineal, comprobar qué fase se duplica y, en consecuencia, retirar o reformular el giro del director en η.

### G7. La energía del nudo: la cota «≲10 eV» no está justificada y los propios datos de positronio la mejoran en dos órdenes de magnitud
**Dónde:** postulado 7, §4, §6 puntos 5 y 6, tabla de calibrados.

- «La ausencia de una torre de estados de rotación exige E_nudo ≲ 10 eV» se afirma sin cálculo. Una estimación de sólido rígido, ΔE ~ (ħc/R)²/E con R = √3 ξ₀ (ħc/R ≈ 0,59 MeV), da ΔE ≈ 35 GeV para E = 10 eV. Estados cargados así estarían al alcance de LEP y del LHC, que buscan electrones excitados hasta escalas de TeV. Habría que derivar la cota y decir con qué datos se compara.
- He rehecho la estimación del positronio: P = (4/3)(√3 ξ₀/a_Ps)³ = 4,2×10⁻⁸ (correcto). Pero el desplazamiento de 0,4 MHz se obtiene con E ≈ 0,04 eV. Con la cota del propio texto (10 eV) saldrían unos **100 MHz**, excluidos por la medida de Ishida et al. (2014). El dato del positronio exige ya E_nudo ≲ 0,1 eV, una cota mucho más fuerte que la de 10 eV, y el texto no lo dice. En el muonio salen 3,2 MHz con 0,04 eV (correcto) y unos 800 MHz con 10 eV.
- Los 0,4 MHz son una **cota superior** sin signo ni valor, no una predicción. Llamarlos «una señal propia del modelo» es excesivo.
- La supresión ≥10³ en el muonio mediante «estrechamiento por movimiento» se introduce para salvar el desacuerdo, y además depende de que la interacción cambie de signo con la orientación, lo que no está justificado. Es un mecanismo *ad hoc*.
- Detalle: dos nudos de radio R se solapan cuando su separación es menor que 2R, no que R. El factor geométrico es 8 veces mayor.

**Propuesta:** sustituir la cota de 10 eV por la que dan el positronio y el muonio. Si se mantiene la propuesta ħc/R_u, decir explícitamente que el muonio ya obliga a un mecanismo de supresión, y etiquetarlo [H].

### G8. Apariencia puntual y *g*−2: el conflicto central se resuelve por postulado
**Dónde:** §3 («Por qué el nudo no tiene carga propia»), §6 puntos 1, 7 y 8, Conclusiones.

El modelo coloca en las extendidas una estructura de radio √3 ξ₀ ≈ 0,87 ƛ_C que porta el eje del espín y en la que está «anclada la fuente del campo». Una subestructura de ese tamaño daría, según Brodsky y Drell, δa ~ m_e R ~ O(1), salvo que haya un mecanismo de protección. El texto lo evita postulando que la fuente «se comporta como puntual en la interacción» [H] y que el nudo no lleva carga. Es exactamente la afirmación que habría que demostrar. Mientras el factor de forma de la circulación no esté calculado, esta es la prueba más fácil de que el modelo falle, y debería presentarse como tal, en el resumen y no solo en §6. Además, la «carga donde está el nudo» del punto 8 (ancla) choca con «toda la carga es circulación» y con que la circulación esté distribuida por la interfaz. Hay que decir dónde está, en las extendidas, la densidad de carga de un electrón aislado.

### G9. El muón como electrón excitado: masa extra, desintegración y contradicciones con el falsador
**Dónde:** Paso 3, §6 punto 11, Conclusiones.

- La composición en cuadratura m² = m_k² + m'² es correcta: un término de masa escalar combinado con el término quiral k γ^η da masa física √(k² + m'²). Pero ese término escalar m' es una **masa en 5D**, y eso contradice el postulado 2 («no pueden tener masa propia en 5D»). Hay que explicar de dónde sale sin propagación en ξ.
- La supresión de μ → eγ se argumenta cualitativamente («difícilmente se convierte»). La cota experimental es de orden 10⁻¹³ en razón de ramificación. Tampoco se discute la conservación del número leptónico (el μ produce siempre ν_μ), que es el dato que una imagen de «electrón excitado» tiene que explicar primero.
- Las conclusiones dan como falsador «un estado del electrón de espín 3/2 o excitado, con su misma carga». El muón es, en el modelo, un estado excitado con su misma carga. Hay que restringirlo («en torno a 3m_ec² y con carga 3e», si se acepta G1).
- El punto 10 dice que «con la carga proporcional a k» el modo 3/2 tendría carga 3e. El Paso 4 y el postulado 5 dicen que la carga la fija la interfaz y es universal. Hay que decidir cuál de las dos vale.

### G10. Etiquetas de estatus que no cumplen las definiciones del propio trabajo
**Dónde:** varias. La nota metodológica define [N] como «reproducible con el código público» y [D] como «se sigue de los postulados».

- **[N] sin código:** tipos II («cálculo no incluido en el repositorio»), masas de los neutrinos (ídem), solapamientos del punto 7 de §2.1 (ídem). Los solapamientos los he reproducido: son |⟨Ψ(η)|Ψ(η−π/2)⟩| = 0,050, 1,6×10⁻⁵ y 1,3×10⁻¹⁵ (no «del orden de 10⁻¹⁶»). Pero el texto no define la medida que usa.
- **[D] apoyado en [H]:** el Paso 2 se etiqueta «[D, orden ½]» cuando el orden ½ es un postulado [H]. El signo de Coulomb es «D por analogía» (véase G5). En la tabla de conclusiones, «*g*_μ ≈ 2: [D]» depende de m' (G9). La «restricción que deja la quiralidad» es [D] dado un [H], y el postulado 4 la presenta como condición impuesta.
- **[O] mal usada:** «Pista especulativa [O]» (punto 12) no es una observación. La coincidencia del 1,6 % con el «trabajo sobre hadrones», que no se ha publicado, no puede comprobarse. El «¾ = s(s+1)» del término de Darwin es una coincidencia numérica y debería ir como [O], no integrada en el argumento.
- **Tercer modo del neutrino:** «masa del orden de las centésimas de eV [H]», cuando «la condición que fija su límite exterior está abierta [A]». Sin esa condición no hay masa.

**Propuesta:** añadir una etiqueta compuesta explícita («D | H₁, H₂», es decir, derivado dadas las hipótesis H₁ y H₂), o rebajar a [H] todo lo que dependa de una [H] no trivial. Retirar [N] de lo que no esté en el repositorio, o subir el código.

### G11. La permeación nula se atribuye a una viscosidad que no puede producirla
**Dónde:** postulado 1, postulado 4, postulado 6, §4 («Por qué el nudo no se desata solo»), §6 punto 2.

El postulado 1 atribuye la permeación nula a una viscosidad **impar** (de Hall), que no disipa y que «desvía» el flujo. Pero una viscosidad impar no impide un flujo a través de las capas: lo rota. Además, el postulado 4 («menos/más viscoso en ξ») y el postulado 6 («la permeación nula y la alta viscosidad en ξ impiden» las inclinaciones) usan lenguaje de viscosidad disipativa. En un esméctico ordinario, la permeación está controlada por un coeficiente disipativo. *g* = 2 depende del postulado 6, y la estabilidad del nudo por analogía con Alfvén depende de la permeación nula, de modo que esta base tiene que ser coherente.

**Propuesta:** formular la permeación nula como una restricción cinemática (una ligadura) y no como un efecto de la viscosidad, y retirar la referencia a la viscosidad impar o limitarla a una posibilidad.

---

## (c) Problemas menores

- **M1.** §6 anuncia «Dieciséis problemas sin resolver», pero la lista tiene **17** puntos (de *g*−2 a «Varias partículas»).
- **M2.** El postulado 1 dice que «η es siempre un ángulo». El Paso 2 dice que «η es la longitud de arco a lo largo de la capa» y usa Ψ(η + 2πξ₀). El punto 1 de la lista vuelve al ángulo con ξ₀⁻¹∂_η. Hay que unificar.
- **M3.** La numeración con `\setcounter{section}{-1}` y `\setcounter{subsection}{-1}` produce un «apartado 2.0» (§2.0) poco habitual. «La relación ξ₀ = ħ/(2mc) del Paso 3» se cita en §2 antes de que se defina el Paso 3 en §2.1.
- **M4.** Las tablas de convergencia y de la región ordenada dan Q ≈ −0,97/−0,98, mientras que el texto habla siempre de Q = 1. Falta una nota sobre la convención de signo.
- **M5.** «Estados de masa, no de sabor» (tabla de §2.0), y a la vez se habla del «modo tauónico» y de ν_τ. Hay que elegir una nomenclatura.
- **M6.** La comparación con los 2,996 ± 0,007 neutrinos de LEP (Janot y Jadach 2020, cifra correcta) se basa en la anchura invisible del Z, que es una interpretación del Modelo Estándar. En un modelo sin Z hay que decir qué se está comparando.
- **M7.** «Un indicio a favor: […] el plano de polarización de una onda que avanza a c (índice 1) no es arrastrado (Player; Jones)». En Player (1976) el arrastre es proporcional a (n_g − 1/n), y se anula trivialmente para n = 1, es decir, sin medio. No es un indicio, y conviene retirarlo.
- **M8.** Bourilkov (2001) existe y da límites al 95 % de 8,2–21,3 TeV (comprobado), pero su resultado principal es un **indicio de 2,6σ** de interacción de contacto axial (Λ ≈ 10,3 TeV). Xi, Liu y Ulbricht (2024) dan desviaciones de más de 3σ en las secciones eficaces diferenciales e indican que «la interacción de contacto podría existir». Resumirlo como «que sus autores no consideran concluyentes» suaviza lo que dicen. Ambos indicios deberían citarse tal cual, porque son relevantes para un modelo con estructura.
- **M9.** Ichinose (2002): su resumen habla de que «aparecen de forma natural» términos de momento dipolar en Kaluza–Klein con formalismo de Cartan. Es aceptable, pero conviene citar la ecuación concreta. Sería útil citar también la literatura de Kaluza–Klein teleparalelo (equivalente teleparalelo de Kaluza–Klein, Andrade, Guillen y Pereira), que es la comparación natural del postulado 6.
- **M10.** Alcerro Mena (2018) está en una revista de baja visibilidad. Se usa como precedente cercano, así que conviene comprobar la referencia (volumen y páginas) y dar un enlace.
- **M11.** Las conclusiones dicen que el nudo «se puede formar desde el vacío por fusión local y un flujo helicoidal», pero §4 aclara que el flujo está impuesto y que el resultado «no [muestra] que el nudo surja espontáneamente del vacío». Hay que alinear las dos frases.
- **M12.** El punto 5(iii) de §6 escribe «κ₂∫c²(ξ₀⁻¹∂_η**n**)²». Sobra c², o hay que definir las dimensiones de κ₂.
- **M13.** La prueba de punto de silla se hace solo con la malla 61³ y caja 5,4, la más gruesa de la tabla. Conviene repetirla en 101³ o 150³.
- **M14.** Las referencias a «la versión 1.0», a «un informe de árbitro» y a lo que «se añadió en la versión 1.1» distraen al lector nuevo (§5, apéndice, nota metodológica). Mejor llevarlas a un registro de cambios separado.
- **M15.** La figura 8 (espín y flujo) no se cita en el apartado «Espín +½ y −½ en el flujo», que es el que la explica. Solo se cita en §6. Las etiquetas de la figura 1 son ilegibles a tamaño de impresión. La figura 9 queda entre la lista de §6 y una frase suelta («Prueba futura…») que parece huérfana.
- **M16.** «En el interior de C\* […] solo caben tres soluciones» se afirma sin demostración, solo con etiqueta [H]. Es una afirmación de completitud y necesita el operador de G3.
- **M17.** El término de Darwin: la equivalencia ⟨δr²⟩ = ¾ƛ_C² = 3ξ₀² y la razón α²/n respecto a la energía de enlace son correctas (comprobado). La frase «Por tanto, esa corrección del orden de α² no es un efecto nuevo…» es redundante: es el término de Darwin por construcción.
- **M18.** En el postulado 1, el tiempo es «un parámetro de evolución externo, […] no una dimensión geométrica». Pero el Paso 1 usa una métrica diag(1,−1,−1,−1,−1) con t dentro del álgebra de Clifford. Es una tensión de presentación que conviene comentar.
- **M19.** En la tabla «Lo que se ha calibrado» faltan parámetros libres que el texto usa: el q del electrón, el umbral de nucleación, K_u y κ, y la profundidad de ablandamiento ε.

**Cifras comprobadas y correctas:** ξ₀ = 1,931×10⁻¹³ m; r_e = 2,818×10⁻¹⁵ m; √3 ξ₀ = 0,866 ƛ_C; u_max/c = 0,0084; ħc/R_u = 0,039 eV (con 5,1 µm); el coste de fundir dos nudos, α ħc/(√3 ξ₀) = 4,3 keV; Koide = 0,666661 (desviación relativa de 9×10⁻⁶, compatible con la incertidumbre de m_τ); (m_τ − m_e)/(m_μ − m_e) = 16,9; a_{1/2} ≈ ¼ − ⅔q²; E/(32√2π²) = 544,35/446,65 = 1,219; 36 % y 18 % de energía liberada en la fusión Q = 2 (con E₂ ≈ 2,00 de Sutcliffe); el álgebra de las matrices γ^A, con (γ⁰γ¹γ²γ³)² = −1; el signo de la transformación gauge; la conversión de Λ ≈ 10 TeV a unos 10⁻²⁰ m; DESI DR2 (Σm_ν < 0,064 y < 0,16 eV); Hanhart et al. (R < 0,66 µm para dos dimensiones); Lee et al. 2020 (52 µm).

**Referencias comprobadas:** Bourilkov 2001, Xi et al. 2024, Hanhart et al. 2001, DESI DR2 e Ichinose 2002 existen y en lo esencial dicen lo que se les atribuye, con los matices de M8 y M9. El resto de referencias clásicas (Frank, Mermin y Ho, Battye y Sutcliffe, Duff, Liu y Rahmfeld, Hosoya et al., Volovik y Mineev, Bressi, Pruttivarasin, Salart, Bancal) son correctas en revista, volumen y contenido según mi conocimiento, aunque no las he revisado todas una por una en la web.

---

## (d) Lo que el trabajo hace bien

- **Honestidad epistémica poco común.** Separa lo calibrado de lo deducido, declara que «no hemos deducido ningún número que no se haya introducido antes como dato» y documenta las vías descartadas y los errores corregidos (§5). Es mejor práctica que la de muchos trabajos profesionales.
- **La parte numérica del hopfión es sólida.** Reproduce la energía de la literatura (1,22) con un análisis de convergencia que distingue el efecto de la caja del de la malla, comprueba el virial (E₂/E₄ → 1) y hace una prueba de punto de silla. Las simulaciones de formación y fusión son originales y se interpretan con cautela («lo robusto es la existencia de un umbral de duración, no su valor»).
- **La reducción a Dirac es correcta en lo técnico.** El papel de la conexión de espín en *g* = 1 frente a *g* = 2 está bien identificado. Aislar el ingrediente no estándar en un único postulado (el 6) es la forma correcta de presentarlo.
- **La lectura de Floquet y Bloch** (circulación en el centro de banda y no en el borde, wronskiano, ausencia de nodos) es correcta y está bien explicada. Los números de la figura 3 son exactos.
- **Ofrece falsadores concretos:** un cuarto neutrino ligero, estados en torno a 3m_ec² y el factor de forma. También discute límites experimentales reales (Eöt-Wash, SN 1987A, LEP, Pruttivarasin, Bressi) en lugar de ignorarlos.
- **Declara el uso de IA** con detalle y responsabilidad.

---

## (e) Valoración global

Como trabajo explicativo para Zenodo es legible, honesto y está mejor documentado de lo habitual en propuestas independientes. Sin embargo, su tesis central (que el espín es una circulación en la dimensión compacta) choca con la propia reducción de Kaluza–Klein que el trabajo usa: allí la circulación compacta es la carga, y el espín ½ viene del álgebra de Clifford (G1, G6). Además, el ingrediente que da la «geometría del espín» (las nucleaciones con |q| grande) anula precisamente la circulación que da la carga y el ħ/2 (G2). Ninguno de los dos problemas es cosmético: afectan al resumen, a las conclusiones y a la exclusión de los «estados de espín 3/2». El sector eléctrico tiene un hueco lógico en el paso de la libertad de origen a «sin término de masa» (G5). Y el sector de neutrinos, que es la única calibración de R_u, no se puede comprobar porque la ecuación del plano compacto no está escrita y el cálculo no está en el repositorio (G3).

Mi recomendación es **revisión mayor antes de publicar**:

1. Reformular la relación espín, carga y ν (G1) y retirar o reformular k_director = 2k_espinor (G6).
2. Calcular ⟨−i∂_η⟩ con el modo completo en función de q (G2).
3. Escribir el operador del plano compacto y publicar el cálculo de los neutrinos (G3).
4. Sustituir la cota de 10 eV por la del positronio (G7).
5. Rebajar las etiquetas que no cumplen las definiciones (G10).

Con esos cambios, el trabajo quedaría como una propuesta especulativa pero coherente. Sin ellos, un físico externo encontrará en la primera lectura la contradicción entre «espín = circulación en η» y «carga = circulación en η».
