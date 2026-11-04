# Guion de Exposición Oral (15 Minutos) y Defensa de Preguntas (5 Minutos)

**Proyecto:** Reducción de artefactos oculares en el EEG de polisomnografía durante el sueño REM mediante descomposición en valores singulares (SVD)  
**Asignatura:** Metodología de la Investigación (Segunda Entrega: Ponencia)  
**Institución:** Universidad Distrital Francisco José de Caldas, Facultad de Ingeniería  
**Docente Asesor:** Prof. Duvan Tellez Castro  
**Autores y Expositores:**
- **Cristian Aguilar** (`20241020155`) : **5 min 00 s** (Diapositivas 1, 2, 7, 8 y 13)
- **Sergio Osorio** (`20241020073`) : **5 min 00 s** (Diapositivas 3, 4, 9, 10 y apertura de 14)
- **Nicolas Velasco** (`20241020009`) : **5 min 00 s** (Diapositivas 5, 6, 11, 12 y cierre de 14)

---

## Guía Rápida de Bolsillo: Cómo explicar las 3 ecuaciones sin enredarse

Si el profesor señala alguna fórmula en la presentación, **no necesitan recitar teoremas extraños**; todo el proyecto se resume en estas **3 ideas sencillas**:

1. **¿Qué es la matriz $X$ y la ecuación $X = U \Sigma V^\top = \sigma_1 u_1 v_1^\top + \sigma_2 u_2 v_2^\top + \sigma_3 u_3 v_3^\top$? (Diapositiva 7)**
   - Cada época de sueño dura $30\text{ segundos}$ a $100\text{ muestras por segundo}$ ($3000\text{ datos}$). Como tenemos $3$ canales (`Fpz-Cz`, `Pz-Oz` y `EOG`), armamos una tabla $X$ de **$3000\text{ filas por }3\text{ columnas}$**.
   - La SVD simplemente separa esa tabla en la suma de **3 componentes independientes (Modo 1, Modo 2 y Modo 3)** ordenados de mayor a menor energía:
     - $u_i$ es **la forma de la onda en el tiempo** ($3000$ puntos).
     - $\sigma_i$ es **qué tanta amplitud o fuerza tiene ese modo** (y su porcentaje de energía es $\sigma_i^2 / (\sigma_1^2 + \sigma_2^2 + \sigma_3^2)$).
     - $v_i$ indica **qué tanto aparece ese modo en cada uno de los 3 canales**.
2. **¿Cómo limpiamos el artefacto ocular? (Diapositiva 7 y 8)**
   - Calculamos la correlación $|r|$ entre la onda temporal de cada modo ($u_i$) y el canal del ojo (`EOG`).
   - Si esa correlación es mayor o igual a $0.50$ ($\tau^* = 0.50$), sabemos que ese modo es el movimiento del ojo, **así que lo apagamos ($\sigma_i = 0$)** y sumamos solo los modos restantes para reconstruir el EEG limpio.
3. **¿Por qué funciona muy bien en microvoltios ($\mu\text{V}$) pero falla con Z-score cuando el artefacto es pequeño ($\alpha = 0.25$)? (Diapositiva 12)**
   - En unidades reales ($\mu\text{V}$), la señal del ojo ($\approx 39.7\,\mu\text{V}$) es casi **4 veces más grande** que el EEG cerebral ($\approx 10.0\,\mu\text{V}$). Como la SVD ordena de mayor a menor tamaño, siempre pone al ojo solo en el **Modo 1** ($\sigma_1 \gg \sigma_2$).
   - Cuando aplicamos **Z-score**, dividimos cada canal por su desviación estándar para que los tres canales tengan el mismo tamaño (varianza $= 1$). Si además el artefacto ocular es pequeño ($\alpha = 0.25$), la onda del ojo y la onda del cerebro quedan empatadas en energía ($\sigma_1 \approx \sigma_2$), por lo que la SVD no logra distinguirlas y las revuelve entre el Modo 1 y el Modo 2.

---

## Mapa General de Tiempos y Transiciones (15:00 Minutos Exactos)

| Bloque de la Guía | Diapositivas | Tiempo Asignado | Minutaje Acumulado | Expositor Responsable |
| :--- | :---: | :---: | :---: | :--- |
| **1. Título y problema de investigación** | Diap. 1 a 3 | **2 min 30 s** | `00:00 - 02:30` | Cristian Aguilar (Diap. 1, 2) $\to$ Sergio Osorio (Diap. 3) |
| **2. Estado del arte y línea temporal** | Diap. 4 a 5 | **3 min 30 s** | `02:30 - 06:00` | Sergio Osorio (Diap. 4) $\to$ Nicolas Velasco (Diap. 5) |
| **3. Prueba de concepto (Datos y SVD)** | Diap. 6 a 7 | **3 min 00 s** | `06:00 - 09:00` | Nicolas Velasco (Diap. 6) $\to$ Cristian Aguilar (Diap. 7) |
| **4. Resultados cuantitativos y cualitativos** | Diap. 8 a 12 | **4 min 00 s** | `09:00 - 13:00` | Cristian (Diap. 8) $\to$ Sergio (Diap. 9, 10) $\to$ Nicolas (Diap. 11, 12) |
| **5. Discusión, conclusiones y código** | Diap. 13 a 14 | **2 min 00 s** | `13:00 - 15:00` | Cristian (Diap. 13) $\to$ Sergio y Nicolas (Diap. 14) |
| **Sesión de Preguntas del Jurado** | Respaldo | **5 min 00 s** | `15:00 - 20:00` | Los 3 integrantes según eje temático |

---

## BLOQUE 1: Título y Problema de Investigación (`00:00 - 02:30` | 2.5 min)

### Diapositiva 1: Portada Institucional (`00:00 - 00:30` | 30 s)
**Expositor: Cristian Aguilar**

> **[Indicaciones escénicas]:** *Mirar al docente y a los compañeros con seguridad. Señalar el título sin leer los códigos numéricos.*

**Parlamento:**  
"Buenos días, profesor Duvan Tellez y compañeros. Nuestro equipo está conformado por Sergio Osorio, Nicolas Velasco y quien les habla, Cristian Aguilar. Hoy presentamos nuestra ponencia titulada *Reducción de artefactos oculares en el EEG de polisomnografía durante el sueño REM mediante descomposición en valores singulares (SVD)*. En estos 15 minutos presentaremos el problema clínico, los ajustes realizados desde el póster, la revisión del estado del arte y nuestra prueba de concepto sobre registros reales de la base de datos *Sleep-EDF*."

---

### Diapositiva 2: Planteamiento del Problema de Investigación (`00:30 - 01:30` | 60 s)
**Expositor: Cristian Aguilar**  
*Título en pantalla: El solapamiento espectral EOG-EEG en banda delta (0.5 a 4 Hz) degrada el EEG frontal en sueño REM.*

> **[Indicaciones escénicas]:** *Señalar el diagrama izquierdo (ojo $\to$ electrodos) y luego el recuadro rojo derecho (dilema del filtro Butterworth).*

**Parlamento:**  
"Durante el sueño REM ocurren al mismo tiempo actividad cerebral rápida, ondas lentas y ráfagas de movimientos oculares rápidos. Como vemos en el diagrama izquierdo, el ojo funciona como un dipolo eléctrico con carga positiva en la córnea y negativa en la retina. Cada vez que el ojo se mueve durante el sueño REM, genera ondas grandes de entre $50$ y $200\text{ microvoltios}$ que viajan por los tejidos de la cabeza y contaminan fuertemente el canal frontal `Fpz-Cz` y, con menor intensidad, el canal posterior `Pz-Oz`.

¿Por qué es tan difícil limpiar este artefacto con un filtro clásico en frecuencia como Butterworth? Porque los movimientos oculares ocurren entre $0.5$ y $4.0\text{ Hz}$, exactamente en la misma banda de frecuencia que las ondas delta cerebrales. Si ponemos un filtro que corte en $0.5\text{ Hz}$, el artefacto ocular pasa completo; y si subimos el corte a $4\text{ Hz}$, borramos la actividad delta real del cerebro. Además, en polisomnografía estándar como *Sleep-EDF*, solo contamos con 2 canales de EEG y 1 canal de EOG, donde métodos como ICA no funcionan por falta de electrodos. Nuestra propuesta es usar la SVD para separar la señal en 3 componentes independientes, aislar el movimiento ocular en el Modo 1 y conservar el EEG limpio en los Modos 2 y 3."

* **[Pase de turno]:** *"A continuación, Sergio explicará cómo aterrizamos este problema desde nuestro póster inicial."*

---

### Diapositiva 3: Evolución desde el Póster, Pregunta, Objetivo e Hipótesis (`01:30 - 02:30` | 60 s)
**Expositor: Sergio Osorio**  
*Título en pantalla: La acotación del póster hacia canales estándar de PSG permite una prueba de concepto verificable.*

> **[Indicaciones escénicas]:** *Señalar brevemente los 4 puntos de la tabla superior y luego las tres tarjetas inferiores.*

**Parlamento:**  
"Gracias, Cristian. En la tabla superior resumimos los cuatro ajustes que hicimos respecto al póster de la primera entrega:
1. Primero, pasamos de querer limpiar múltiples artefactos a enfocarnos **únicamente en el artefacto ocular durante el sueño REM**, donde el choque en la banda delta es crítico.
2. Segundo, en lugar de suponer cascos de 19 o más electrodos, aterrizamos el método a los **3 canales reales de *Sleep-EDF*** (`Fpz-Cz`, `Pz-Oz` y `EOG horizontal`), armando matrices de $3000\text{ muestras por }3\text{ columnas}$ en cada época de 30 segundos.
3. Tercero, como en un EEG clínico nunca conocemos cuál era la onda cerebral exacta antes de que se ensuciara, diseñamos una **mezcla controlada con datos reales**: tomamos épocas REM tranquilas libres de parpadeos como referencia limpia ($x_{\text{ref}}$) y les sumamos movimientos oculares reales en tres niveles de intensidad ($\alpha = 0.25, 0.50\text{ y }0.80$).
4. Y cuarto, definimos una regla automática basada en correlación ($|r| \ge 0.50$) calibrada en Desarrollo y comparamos qué pasa al usar los datos en microvoltios frente a normalizarlos con Z-score.

De ahí surgen nuestra pregunta, nuestro objetivo general y nuestra hipótesis: que la SVD logra una ganancia de relación señal a ruido mayor que $0\text{ dB}$, con un intervalo de confianza al $95\%$ positivo en los sujetos de evaluación."

---

## BLOQUE 2: Estado del Arte y Línea Temporal (`02:30 - 06:00` | 3.5 min)

### Diapositiva 4: Estado del Arte I - Diagrama de Relaciones (`02:30 - 04:15` | 1 min 45 s)
**Expositor: Sergio Osorio**  
*Título en pantalla: La literatura evoluciona desde regresión lineal e ICA hacia filtrado por componentes con SVD.*

> **[Indicaciones escénicas]:** *Señalar los 4 bloques superiores de izquierda a derecha y cómo bajan hacia el bloque azul de nuestro trabajo.*

**Parlamento:**  
"Para ubicar nuestro trabajo en la literatura, organizamos las 13 referencias en cuatro grupos:
- En el **primer grupo**, Berg y Picton `[3]` y Gratton *et al.* `[2]` propusieron restar del EEG una fracción fija del canal EOG mediante regresión lineal, mientras Kemp *et al.* `[1]` publicaron la base de datos *Sleep-EDF*. El problema de la regresión simple es que el electrodo ocular también capta algo de actividad cerebral frontal, por lo que al restarlo directamente se deforma el EEG.
- En el **segundo grupo**, Makeig `[4]`, Vigário `[5]` y Jung `[6]` introdujeron la separación ciega de fuentes mediante ICA. ICA funciona muy bien cuando tenemos 16, 32 o 64 electrodos, pero falla cuando solo tenemos 2 canales de EEG y 1 de EOG.
- En el **tercer grupo**, Golub y Van Loan `[7]`, Lagerlund `[8]` y Gavish y Donoho `[9]` fundamentaron el uso de la descomposición en valores singulares (SVD). Lagerlund demostró en 1997 un punto clave para nuestro proyecto: la SVD separa muy bien un artefacto ocular siempre que ese artefacto tenga mucha más amplitud que el EEG de fondo.
- Y en el **cuarto grupo**, autores recientes entre 2021 y 2025 como Salaman `[11]`, Grochala `[12]`, Ghosh `[10]` y Phadikar `[13]` combinan SVD con transformadas Wavelet, pero requieren ajustar muchos parámetros y consumen mucho más tiempo de cómputo."

* **[Pase de turno]:** *"Ahora Nicolas presentará la línea de tiempo y el vacío que resuelve nuestra prueba de concepto."*

---

### Diapositiva 5: Estado del Arte II - Línea Temporal y Vacío Identificado (`04:15 - 06:00` | 1 min 45 s)
**Expositor: Nicolas Velasco**  
*Título en pantalla: Existe un vacío metodológico en PSG de pocos canales entre el filtrado clásico y los métodos complejos.*

> **[Indicaciones escénicas]:** *Recorrer la línea temporal de izquierda (1978) a derecha (2026) y luego señalar las dos tarjetas inferiores.*

**Parlamento:**  
"Gracias, Sergio. En esta línea de tiempo vemos el recorrido desde 1978 hasta llegar a nuestro trabajo en 2026: desde la regresión de Gratton `[2]` en 1983, pasando por el primer filtro SVD de Lagerlund `[8]` en 1997 y la base *Sleep-EDF* de Kemp `[1]` en el año 2000, hasta los métodos híbridos de Phadikar `[13]` en 2025.

Al revisar estos antecedentes encontramos **tres vacíos claros**:
1. Primero, los métodos más famosos como ICA `[4]`, `[5]` no se pueden aplicar en estudios básicos de sueño que solo tienen 2 canales de EEG.
2. Segundo, casi ningún trabajo mide al mismo tiempo cuánto mejora la relación señal a ruido y qué tan bien se recupera la banda delta ($0.5\text{ a }4\text{ Hz}$) durante el sueño REM.
3. Y tercero, no se había explicado de forma práctica qué pasa si antes de aplicar la SVD cometemos el error de normalizar los canales con Z-score en vez de dejarlos en microvoltios.

Nuestro aporte consiste en un método con SVD directo, que tarda menos de un milisegundo por época usando solo 3 canales, valida sus resultados en pacientes nuevos y demuestra por qué mantener la escala en microvoltios es indispensable."

---

## BLOQUE 3: Prueba de Concepto (`06:00 - 09:00` | 3.0 min)

### Diapositiva 6: Prueba de Concepto I - Flujo Metodológico y Datos Reales (`06:00 - 07:30` | 1 min 30 s)
**Expositor: Nicolas Velasco**  
*Título en pantalla: El diseño semi-sintético sobre 10 sujetos de Sleep-EDF garantiza verdad de referencia libre de fugas.*

> **[Indicaciones escénicas]:** *Señalar los 5 pasos de la Figura 1 a la izquierda y la tabla de los 10 pacientes a la derecha.*

**Parlamento:**  
"Para poner a prueba el método, implementamos el flujo de cinco etapas de la Figura 1 usando datos reales de *Sleep-EDF Expanded v1.0.0* `[1]`, grabados a $100\text{ Hz}$.

Descargamos todas las épocas de sueño REM de **10 pacientes completos**, sumando **1,594 épocas de 30 segundos** (más de 13 horas de sueño REM). Para hacer un experimento limpio y sin trampa, dividimos los 10 pacientes en dos grupos totalmente separados, como muestra la tabla derecha:
- **Grupo de Desarrollo (`SC4001` a `SC4041`, 5 pacientes):** 863 épocas REM, de donde seleccionamos $N_{\text{dev}} = 282$ épocas limpias sin movimientos oculares ($|r| < 0.15$). Este grupo se usó únicamente para calibrar el umbral del algoritmo.
- **Grupo de Evaluación (`SC4051` a `SC4091`, 5 pacientes nuevos):** 731 épocas REM, de donde obtuvimos $N_{\text{eval}} = 183$ épocas limpias para la prueba final.

A cada época limpia ($x_{\text{ref}}$) le sumamos ráfagas reales de movimiento ocular del mismo paciente: sumamos $\alpha \cdot eog(t)$ en el canal frontal `Fpz-Cz` y el $30\%$ de ese valor ($0.30\,\alpha \cdot eog(t)$) en el canal posterior `Pz-Oz`, porque está más lejos de los ojos. Probamos tres niveles de contaminación: baja ($\alpha=0.25$), moderada ($\alpha=0.50$) y severa ($\alpha=0.80$)."

* **[Pase de turno]:** *"Le doy paso a Cristian para explicar de forma sencilla cómo opera la SVD sobre cada época."*

---

### Diapositiva 7: Prueba de Concepto II - Formulación Matricial SVD y Línea Base (`07:30 - 09:00` | 1 min 30 s)
**Expositor: Cristian Aguilar**  
*Título en pantalla: La SVD descompone la época de 3 canales en 3 modos independientes y apaga el modo ocular.*

> **[Indicaciones escénicas]:** *Apuntar a la suma de los 3 modos en la tarjeta superior izquierda y luego a los dos regímenes de preprocesamiento a la derecha.*

**Parlamento:**  
"Gracias, Nicolas. ¿Cómo funciona exactamente nuestro filtro SVD en cada época de 30 segundos? Como la señal tiene $100\text{ muestras por segundo}$, en 30 segundos tenemos $3000\text{ datos}$ por cada canal. Ponemos los 3 canales (`Fpz-Cz`, `Pz-Oz` y `EOG`) uno al lado del otro formando una tabla $X$ de **$3000\text{ filas por }3\text{ columnas}$**.

Al aplicar la SVD, esa tabla se separa en la suma de **3 componentes o modos independientes** ($\sigma_1 u_1 v_1^\top + \sigma_2 u_2 v_2^\top + \sigma_3 u_3 v_3^\top$), donde cada pieza tiene un significado muy simple:
- $u_i$ es **la forma de la onda en el tiempo** ($3000$ puntos).
- $\sigma_i$ es **la amplitud o fuerza de ese modo**, ordenados siempre del más grande ($\sigma_1$) al más pequeño ($\sigma_3$).
- Y $v_i$ nos dice **cuánto aparece ese modo en cada uno de los 3 canales**.

Entonces, la regla para limpiar el artefacto es muy directa: medimos la correlación $|r|$ de cada uno de los 3 modos con el canal ocular `EOG`. Si la correlación es mayor o igual a $0.50$, sabemos que ese modo es el movimiento ocular, **lo apagamos haciendo su peso $\sigma_i = 0$**, y sumamos los modos restantes para recuperar el EEG limpio.

Para ver el efecto del preprocesamiento, probamos dos variantes:
1. **`SVD_Fisico_uV`**, donde solo restamos el promedio y dejamos las señales en microvoltios reales (donde el ojo mide unos $39.7\,\mu\text{V}$ y el EEG unos $10\,\mu\text{V}$).
2. Y **`SVD_Zscore`**, donde dividimos cada canal por su desviación estándar para que los tres tengan varianza igual a 1.

Comparamos ambos contra la señal sin procesar y contra un filtro **Butterworth pasa-banda de $0.5\text{ a }30\text{ Hz}$**, midiendo la relación señal a ruido ($\text{SNR}$ en decibelios, donde más alto es mejor) y los errores de reconstrucción de la onda y de la banda delta (donde más cercano a cero es mejor)."

---

## BLOQUE 4: Resultados Cuantitativos y Cualitativos (`09:00 - 13:00` | 4.0 min)

### Diapositiva 8: Resultados Cuantitativos I - Calibración del Umbral en Desarrollo (`09:00 - 10:00` | 60 s)
**Expositor: Cristian Aguilar**  
*Título en pantalla: El barrido en desarrollo identifica un umbral óptimo $\tau^* = 0.50$ con una meseta estable en escala física.*

> **[Indicaciones escénicas]:** *Señalar la línea recta verde en el panel derecho de la Figura 3 y la fila $\tau^* = 0.50$ en la tabla.*

**Parlamento:**  
"En la Figura 3 vemos el primer resultado cuantitativo: cómo elegimos el umbral de correlación $\tau^*$ usando únicamente las $282$ épocas de los 5 pacientes de Desarrollo ($\alpha=0.50$).

Si miran el panel derecho de la Figura 3 y la tabla, cuando trabajamos en **microvoltios reales** (línea verde), el método es sumamente estable: da exactamente la misma ganancia de **$+10.89\text{ dB}$** sin importar si ponemos el umbral en $0.30, 0.40, 0.50, 0.60\text{ o }0.70$, apagando siempre $1.01$ modos en promedio (es decir, apaga únicamente el Modo 1). ¿Por qué? Porque el Modo 1 tiene correlación casi de $1.0$ con el ojo y el Modo 2 tiene correlación menor a $0.16$, así que cualquier umbral intermedio los separa sin problema.

En cambio, en el panel izquierdo (con **Z-score**), la curva sube y baja fuertemente y alcanza su mejor valor ($+5.90\text{ dB}$) justo en $\tau^* = 0.50$. Por eso dejamos fijo $\tau^* = 0.50$ para todas las pruebas siguientes."

* **[Pase de turno]:** *"Ahora Sergio presentará los resultados obtenidos en los 5 pacientes nuevos de Evaluación."*

---

### Diapositiva 9: Resultados Cuantitativos II - Energía y Correlación en Evaluación (`10:00 - 11:00` | 60 s)
**Expositor: Sergio Osorio**  
*Título en pantalla: El primer modo singular concentra entre 82.3% y 90.3% de la energía con correlación EOG de 0.98.*

> **[Indicaciones escénicas]:** *Señalar las barras azules altas de la Figura 2 y las tres primeras filas de la tabla derecha.*

**Parlamento:**  
"Gracias, Cristian. Al pasar a los 5 pacientes nuevos de Evaluación ($183$ épocas REM), la Figura 2 y la tabla muestran qué porcentaje de energía queda en cada uno de los 3 modos:
- Cuando usamos **escala física en microvoltios** (barras azules y tres primeras filas de la tabla), el **Modo 1** concentra el **$82.3\%$** de la energía en contaminación baja, el **$86.0\%$** en moderada y el **$90.3\%$** en severa. Además, la correlación del Modo 1 con el canal ocular es de **$0.97\text{ a }0.98$**, mientras que los Modos 2 y 3 tienen correlaciones menores o iguales a **$0.17$**. Esto confirma que todo el artefacto ocular quedó atrapado en el Modo 1.
- En cambio, miren la cuarta fila de la tabla (**Z-score con contaminación baja $\alpha=0.25$**): la energía se parte entre los tres modos ($48.7\%$, $30.9\%$ y $20.4\%$) y la correlación con el ojo se riega en los tres componentes ($0.66$ y $0.41$)."

---

### Diapositiva 10: Resultados Cuantitativos III - SNR, Errores e Intervalos de Confianza (`11:00 - 12:00` | 60 s)
**Expositor: Sergio Osorio**  
*Título en pantalla: La SVD en escala física logra ganancias de +13.9 a +14.9 dB y verifica la hipótesis con IC 95% > 0.*

> **[Indicaciones escénicas]:** *Señalar las barras altas azules en el panel izquierdo de la Figura 4 y los intervalos a la derecha de la línea punteada roja en el panel derecho.*

**Parlamento:**  
"En la Figura 4 y en la tabla derecha vemos el desempeño final sobre el canal frontal `Fpz-Cz` en Evaluación:
1. Primero, el filtro **Butterworth ($0.5\text{ a }30\text{ Hz}$)** casi no limpia el artefacto: en contaminación moderada ($\alpha=0.50$), la $\text{SNR}$ apenas sube de $-2.55\text{ dB}$ a $+1.07\text{ dB}$, dejando un error enorme en la banda delta de $1.263$.
2. Segundo, nuestra **SVD en microvoltios** sube la $\text{SNR}$ hasta **$+18.37\text{ dB}$** en contaminación baja, **$+11.39\text{ dB}$** en moderada y **$+7.35\text{ dB}$** en severa. En $\alpha=0.50$, baja el error de la onda (Frobenius) de $1.305$ a **$0.257$** y reduce el error en la banda delta de $1.263$ a **$0.419$** (un $66.8\%$ menos de error que Butterworth).
3. Y tercero, en el panel derecho de la Figura 4 vemos los intervalos de confianza al $95\%$ calculados por paciente con Bootstrap: para la **SVD en microvoltios**, todos los intervalos quedan muy a la derecha de la línea roja de $0\text{ dB}$ (ganando entre $+9.1\text{ y }+14.5\text{ dB}$ sobre Butterworth), **comprobando completamente nuestra hipótesis**. En cambio, el Z-score en contaminación baja ($\alpha=0.25$) cruza la línea roja ($[-1.63, +2.07]\text{ dB}$), por lo que no logra mejora significativa en ese caso."

* **[Pase de turno]:** *"Le doy la palabra a Nicolas para ver las señales en el tiempo y entender de forma sencilla por qué falla el Z-score."*

---

### Diapositiva 11: Resultados Cualitativos I y II - Modos y Reconstrucción Temporal (`12:00 - 12:30` | 30 s)
**Expositor: Nicolas Velasco**  
*Título en pantalla: El modo dominante captura la morfología ocular en baja frecuencia y restaura la banda delta frontal.*

> **[Indicaciones escénicas]:** *Señalar el Modo 1 rojo en la Figura 5 (izquierda) y la superposición de la línea verde con la negra en la Figura 6 (derecha).*

**Parlamento:**  
"Gracias, Sergio. En la Figura 5 (a la izquierda) vemos visualmente los 3 modos en una época real del paciente `SC4071`: el **Modo 1** (en rojo) captura exactamente las sacadas lentas del ojo ($98.8\%$ de la energía y correlación de $1.00$), mientras que los **Modos 2 y 3** (en azul y verde) muestran únicamente las oscilaciones rápidas del cerebro libres de parpadeos.

A la derecha, en la Figura 6, vemos la reconstrucción de la onda: entre los segundos 14 y 21, cuando ocurren movimientos oculares fuertes, el filtro Butterworth (línea naranja) se deforma igual que la señal sucia. En cambio, la señal reconstruida con SVD (línea verde, $+12.9\text{ dB}$) sigue casi perfecta a la señal cerebral limpia de referencia (línea negra), y en la gráfica de frecuencias quita el exceso falso entre $0.5\text{ y }4\text{ Hz}$ sin tocar el resto del espectro."

---

### Diapositiva 12: Resultados Cualitativos III - ¿Cuándo Funciona y Cuándo Falla? (`12:30 - 13:00` | 30 s)
**Expositor: Nicolas Velasco**  
*Título en pantalla: La normalización Z-score falla en baja contaminación ($\alpha = 0.25$) porque iguala el tamaño de las señales.*

> **[Indicaciones escénicas]:** *Señalar el Panel A (éxito) y Panel B (fallo) en la Figura 7 y luego los 3 pasos numerados de la tarjeta roja derecha.*

**Parlamento:**  
"En la Figura 7 explicamos de forma muy sencilla **cuándo funciona bien nuestro método y por qué falla el Z-score en contaminación baja ($\alpha=0.25$)**, siguiendo los tres puntos del recuadro derecho:
1. **Primero**, cuando dejamos los datos en **microvoltios reales**, el canal del ojo (`EOG`, $\approx 39.7\,\mu\text{V}$) es casi **4 veces más grande** que el canal del cerebro (`EEG`, $\approx 10.0\,\mu\text{V}$). Como la SVD ordena siempre de mayor a menor tamaño, el ojo queda aislado en el Modo 1 ($\sigma_1$ mucho mayor que $\sigma_2$).
2. **Segundo**, cuando aplicamos **Z-score**, dividimos cada canal por su desviación estándar, lo que encoge el canal del ojo y agranda el del cerebro hasta obligarlos a tener el mismo tamaño (varianza igual a 1).
3. **Tercero**, si además la contaminación ocular es pequeña ($\alpha=0.25$, Panel B), la onda del ojo y la onda del cerebro quedan empatadas en tamaño ($\sigma_1 \approx \sigma_2$). Al tener la misma energía, la SVD ya no puede distinguirlas y mezcla el ojo entre el Modo 1 y el Modo 2; cuando el algoritmo apaga ambos modos, termina borrando parte de la señal cerebral en el $47.5\%$ de las épocas."

* **[Pase de turno]:** *"Para cerrar, Cristian presentará la discusión, las limitaciones y el trabajo futuro."*

---

## BLOQUE 5: Discusión, Conclusiones y Código Reproducible (`13:00 - 15:00` | 2.0 min)

### Diapositiva 13: Discusión, Limitaciones y Trabajo Futuro (`13:00 - 14:00` | 60 s)
**Expositor: Cristian Aguilar**  
*Título en pantalla: La SVD supera el filtrado lineal sin entrenamiento, pero requiere ventanas adaptativas.*

> **[Indicaciones escénicas]:** *Recorrer las tres columnas de izquierda a derecha.*

**Parlamento:**  
"Gracias, Nicolas. En la discusión destacamos que, frente al filtro Butterworth que no puede separar dos señales de la misma frecuencia y frente a ICA que exige muchos electrodos `[5]`, la SVD limpia eficazmente el artefacto ocular usando solo los 3 canales básicos de *Sleep-EDF* `[1]`, confirmando la regla práctica de Lagerlund `[8]`: para que la SVD separe el artefacto, hay que conservar su mayor amplitud natural en microvoltios.

De forma autocrítica, señalamos **tres limitaciones de nuestro estudio**:
1. Primero, nuestra mezcla suma el artefacto de forma instantánea con un factor fijo $\alpha$, sin simular posibles retardos a través del cráneo.
2. Segundo, aplicamos una sola SVD a toda la ventana de 30 segundos, aunque las sacadas oculares duren solo unos pocos segundos dentro de la época.
3. Y tercero, nuestra regla necesita que el electrodo ocular `EOG` esté bien conectado.

Como **trabajo futuro**, proponemos aplicar la SVD en **ventanas cortas deslizantes de 2 a 5 segundos** solo cuando aparezca una ráfaga ocular; identificar el modo ocular mirando cuánta energía lenta ($0.5\text{ a }4\text{ Hz}$) concentra en caso de que el electrodo EOG se desconecte; y probar cuánto mejora la clasificación automática de etapas de sueño en toda la base *Sleep-EDF*."

---

### Diapositiva 14: Conclusiones, Referencias IEEE y Repositorio Reproducible (`14:00 - 15:00` | 60 s)
**Expositores: Sergio Osorio (`14:00 - 14:30`) y Nicolas Velasco (`14:30 - 15:00`)**  
*Título en pantalla: La SVD en escala física es una solución eficaz y reproducible para PSG de pocos canales.*

> **[Indicaciones escénicas]:** *Sergio presenta las conclusiones 1 y 2; Nicolas presenta la conclusión 3, las referencias IEEE y el enlace al cuaderno reproducible.*

**Parlamento (Sergio Osorio, `14:00 - 14:30`):**  
"Para concluir:
1. **Primero**, comprobamos nuestra hipótesis en los 5 pacientes independientes de Evaluación ($183$ épocas REM): la SVD en microvoltios con $\tau^*=0.50$ alcanzó relaciones señal a ruido de $+18.37\text{ dB}$, $+11.39\text{ dB}$ y $+7.35\text{ dB}$, superando a la señal sin procesar en más de $+13.9\text{ dB}$ y al filtro Butterworth en $+9.1\text{ a }+14.5\text{ dB}$ con un $95\%$ de confianza.
2. **Segundo**, logramos limpiar tanto el canal frontal `Fpz-Cz` como el posterior `Pz-Oz` recuperando la banda delta cerebral y reduciendo en un $66.8\%$ el error espectral frente a Butterworth."

**Parlamento (Nicolas Velasco, `14:30 - 15:00`):**  
"3. Y **tercero**, demostramos que en montajes de pocos canales con referencia ocular **no se debe normalizar con Z-score**, sino conservar la escala en microvoltios para que el movimiento ocular domine el Modo 1 y no se mezcle con el EEG.

En el recuadro inferior izquierdo dejamos el enlace a nuestro repositorio y al cuaderno de Jupyter (`prueba_concepto_svd_rem.ipynb`), que descarga los datos de PhysioNet y reproduce con semilla fija (`seed = 2026`) todas las tablas y las 7 figuras de esta presentación. Muchas gracias por su atención; quedamos atentos a sus preguntas."

---

## BANCO DE DEFENSA PARA LA SESIÓN DE PREGUNTAS (`15:00 - 20:00` | 5 Minutos)

Respuestas claras, directas y sin matemáticas complicadas ante las preguntas típicas del profesor:

### 1. *"¿Por qué hicieron una mezcla sumando EEG limpio + EOG en vez de aplicar la SVD directamente a una época sucia?"*
- **Respuesta sencilla (Sergio Osorio):**  
  "Porque si tomamos una época que ya viene sucia del hospital, no sabemos cuánto valía la onda cerebral limpia antes del parpadeo. Sin conocer la señal limpia original ($x_{\text{ref}}$), es imposible calcular la fórmula de $\text{SNR}$ o saber si el filtro borró también actividad cerebral. En cambio, al tomar épocas reales de sueño REM donde el ojo estuvo quieto ($x_{\text{ref}}$) y sumarles movimientos oculares reales del mismo paciente, trabajamos con ondas 100% reales pero conocemos la respuesta exacta para medir en decibelios qué tan bien funcionó el filtro."

### 2. *"Explíqueme con sus palabras qué son $U$, $\Sigma$ y $V^\top$ en la ecuación $X = U \Sigma V^\top$."*
- **Respuesta sencilla (Cristian Aguilar):**  
  "Claro, profe. Nuestra matriz $X$ tiene $3000\text{ filas}$ (los puntos en el tiempo de los 30 segundos) y $3\text{ columnas}$ (los 3 canales: frontal, posterior y ojo). La SVD desarma esa tabla en 3 modos independientes:
  - $U$ tiene las **formas de onda en el tiempo** ($3000$ puntos cada una): el Modo 1 tiene la forma del movimiento ocular lento, y los Modos 2 y 3 tienen las ondas rápidas del cerebro.
  - $\Sigma$ (Sigma) es una diagonal con 3 números ($\sigma_1, \sigma_2, \sigma_3$) que nos dicen **qué tan grande o qué tanta energía tiene cada modo**, ordenados de mayor a menor.
  - Y $V^\top$ es una tabla pequeña de $3\times 3$ que nos dice **en qué proporción se reparte cada modo entre los 3 electrodos**. Cuando hacemos $\sigma_1 = 0$, apagamos el Modo 1 del ojo y al multiplicar de vuelta solo quedan los Modos 2 y 3 del cerebro."

### 3. *"¿Por qué separaron 5 pacientes para Desarrollo y 5 pacientes para Evaluación?"*
- **Respuesta sencilla (Nicolas Velasco):**  
  "Para no hacer trampa ajustando el umbral con los mismos pacientes donde íbamos a medir el resultado final. Usamos los primeros 5 pacientes (`SC4001` a `SC4041`, $282$ épocas) para probar distintos valores de umbral entre $0.15\text{ y }0.85$; ahí vimos que $\tau^* = 0.50$ era el punto óptimo. Congelamos ese $0.50$ y luego lo probamos a ciegas en los otros 5 pacientes nuevos (`SC4051` a `SC4091`, $183$ épocas) para demostrar que el método sí generaliza a personas distintas."

### 4. *"¿Qué significa el intervalo de confianza Bootstrap al 95% que muestran en la Figura 4?"*
- **Respuesta sencilla (Sergio Osorio):**  
  "Es una prueba estadística para estar seguros de que la ganancia de decibelios no fue casualidad de unas pocas épocas. El método Bootstrap toma nuestros pacientes de evaluación, los sortea con reemplazo $10\,000$ veces en el computador y calcula en qué rango cae la mejora promedio el $95\%$ de las veces. Como para la SVD en microvoltios todo ese rango cae entre $+7.9\text{ y }+11.3\text{ dB}$ por encima de Butterworth (muy por encima de cero), queda comprobado estadísticamente que nuestro método sí mejora la señal."
