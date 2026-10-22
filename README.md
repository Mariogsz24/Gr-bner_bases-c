# Grobner_bases-c
Algoritmo de Buchberger para obtener una base de Gröbner reducida y minimal de un conjunto de polinomios multivariables.
## Estructura del Proyecto

**Grobner.c**: Implementación principal con asignación de memoria dinámica en el heap para la gestión de listas polinomiales.

**Grobnergmp.c**: Implementación con la libreria GMP para asegurarnos no perder precisión al truncar con un epsilon.

## Características Técnicas

**Órdenes Monomiales Soportados:** Lexicográfico (LEX), Lexicográfico Graduado (GRLEX) y Lexicográfico Graduado Inverso (GREVLEX).
**Salida:** Renderizado automático de polinomios en consola utilizando variables alfabéticas estándar.
**Aritmética:** Procesamiento numérico basado en coma flotante de doble precisión (`double`) con control de tolerancia (`EPSILON`).


## Compilación y Ejecución

Se requiere un entorno de desarrollo con el compilador GCC (Linux, WSL o MSYS2).
Para compilar y ejecutar :

```bash
gcc Grobner.c -o grobner

./grobner

```

Se requiere sustituir el tipo `double` por aritmética de precisión múltiple (ej. librería GMP) para evitar la cancelación de coeficientes


La implementación normal es una prueba de concepto. Para entornos de producción o ideales de alta complejidad, se requiere sustituir el tipo `double` por aritmética de precisión múltiple (ej. librería GMP) para evitar la cancelación de coeficientes

```bash

sudo apt update && sudo apt install build-essential libgmp-dev

gcc Grobner.c -o grobner -lgmp

./grobner

```

