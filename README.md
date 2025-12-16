# Proyecto Multi-CountSketch

Este proyecto implementa una herramienta eficiente en C++ para el análisis de frecuencias de K-mers en secuencias de ADN (formato FASTA) utilizando estructuras probabilísticas (CountSketch). Incluye scripts de visualización en Python para analizar la distribución de frecuencias y detectar anomalías en cromosomas.

* **Integrantes:**
  * Nicolás López
  * Benjamín Espinoza
* **Asignatura:** Topicos en Manejo de Grandes Volumenes de Datos
* **Semestre:** 2025-2
* **Profesora:** Cecilia Hernández

---

## Instalación y Requisitos

### Requisitos previos
* **Compilador C++:** Compatible con C++17 o superior.
* **Python 3.x:** Para la generación de gráficos.
* **Librerías de Python:** `pandas`, `matplotlib`.

### Compilación
Para compilar el proyecto, desde el directorio raíz, ejecutar el siguiente comando:

```bash
g++ -O3 -fopenmp -march=native main_genome.cpp -o mcsketch
```

---

## Uso de la Herramienta

El programa principal `mcsketch` funciona mediante línea de comandos y tiene tres modos distintos.

### Sintaxis General
```bash
./mcsketch <modo> -k <lista_k> -d <dimension> -w <hashes> [-p <pesos>]
```

### Argumentos
* `<modo>`: 
  * `count`: Solo procesa archivos y guarda la estructura (`.bin`).
  * `score`: Carga una estructura existente y calcula puntajes.
  * `both`: Entrena y calcula puntajes en una sola ejecución.
* `-k`: Lista de longitudes de K-mers separadas por comas (ej: `15,21,31`).
* `-d`: Dimensión del Sketch (columnas). Para genoma humano se recomienda 67108864 (2^26).
* `-w`: Ancho/Profundidad del Sketch (filas/hashes). Recomendado: 5.
* `-p` (Opcional): Pesos para el scoring (ej: `1.0,1.0,2.0`), debe tener la misma dimensión que K, sigue el mismo orden.

### Ejemplos de Ejecución

**1. Ejecución completa (Conteo + Scoring):**
Procesa los archivos usando K-mers de longitud 15, 21 y 31.
```bash
./mcsketch both -k 15,21,31 -d 67108864 -w 5
```

**2. Solo calcular puntajes (con pesos personalizados):**
Si ya existe el archivo binario, puedes recalcular puntajes dando más peso a los K-mers más largos.
```bash
./mcsketch score -k 15,21,31 -d 67108864 -w 5 -p 1.0,1.0,2.0
```

---

## Visualización de Resultados

El proyecto incluye un script en Python para generar gráficos comparativos de los resultados obtenidos.

Para generar los gráficos, ejecute desde la raíz:

```bash
python plots/grapher.py
```

Esto generará una imagen en la carpeta `plots/results/`:
* `resultados_score.png`: Histograma con los puntajes por cromosoma.
---

## Set de Datos de Prueba

El repositorio incluye un set de datos de prueba ubicado en la carpeta `datasets/`.

* Estos archivos corresponden a cromosomas completos del genoma humano (GRCh38) en formato FASTA.
* El programa detectará automáticamente todos los archivos `.fa` o `.fasta` en esta carpeta para su procesamiento.