#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <random> // Para generar semillas de hash

// Definimos el tipo de contador
using CounterType = uint64_t; // Soporta frecuencias de hasta ~4 mil millones

class CountSketch {
private:
    const int W;
    const int D;
    std::vector<std::vector<CounterType>> matrix; 

    // Semillas para las funciones de hash. Una para cada fila 'w' para independencia.
    std::vector<uint64_t> seeds_h; // Para h(x) -> columna
    std::vector<uint64_t> seeds_g; // Para g(x) -> signo (+1 o -1)

    /**
 * @brief Función de Hash rápida (MurmurHash3 Finalizer adaptado).
 * * Esta es una función de mezcla que toma una clave de 64 bits (el k-mer) 
 * y aplica una serie de multiplicaciones y XOR-shifts para producir un 
 * hash de 64 bits bien distribuido.
 *
 * @param kmer La clave de 64 bits (el k-mer codificado).
 * @param seed La semilla de 64 bits (para independencia).
 * @return El valor de hash de 64 bits.
 */
uint64_t fast_hash(uint64_t kmer, uint64_t seed) const {
    // Constantes de MurmurHash3 para la mezcla de 64 bits
    const uint64_t C1 = 0x87c37b91114253d5ULL;
    const uint64_t C2 = 0x4cf5ad432745937fULL;
    
    uint64_t h = kmer ^ seed; // Inicializar con la clave y la semilla
    
    // --- Etapa de Mezcla (similar al finalizer de MurmurHash3) ---
    
    // Mezclar con C1 y rotación (similar a la mezcla del bloque de datos)
    h ^= h >> 27;
    h *= C1;
    h ^= h >> 27;
    h *= C2;
    
    // Mezcla final (avalancha) para asegurar una buena distribución
    h ^= h >> 33; 
    h *= 0xff51afd7ed558ccdULL; 
    h ^= h >> 33; 
    h *= 0xc4ceb9fe1a85ec53ULL; 
    h ^= h >> 33;

    return h;
}

public:
    /**
     * @brief Constructor de CountSketch.
     * @param w Profundidad (filas).
     * @param d Ancho (columnas).
     */
    CountSketch(int w, int d) : W(w), D(d) {
        // Inicializar la matriz con W filas y D columnas, todas a 0
        matrix.resize(W, std::vector<CounterType>(D, 0));
        
        // Inicializar las semillas de hash para cada fila
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> distrib;

        for (int i = 0; i < W; ++i) {
            seeds_h.push_back(distrib(gen));
            seeds_g.push_back(distrib(gen));
        }
    }

    /**
     * @brief Incrementa el contador para un k-mer dado.
     * NOTA: En la fase de paralelización, esta función debe ser llamada
     * con mecanismos de bloqueo (locks) si la matriz es compartida.
     * @param kmer El k-mer codificado (uint64_t).
     */
    void update(uint64_t kmer) {
        for (int i = 0; i < W; ++i) {
            // 1. Obtener el índice de columna (0 a D-1)
            uint64_t hash_h = fast_hash(kmer, seeds_h[i]);
            int column_index = hash_h % D;

            // 2. Obtener el signo (+1 o -1)
            uint64_t hash_g = fast_hash(kmer, seeds_g[i]);
            int sign = (hash_g & 1) ? 1 : -1; // Usa el bit menos significativo

            // 3. Actualizar el contador. 
            // En un Count Sketch estándar, se agrega 'sign' al contador. 
            // Para conteo de frecuencias (siempre +1): simplemente incrementamos.
            // Usaremos el incremento simple por ser un 'Count Sketch' aplicado a frecuencias de k-mer.
            matrix[i][column_index]++;
        }
    }

    /**
     * @brief Estima la frecuencia de un k-mer.
     * La estimación es la mediana de las W entradas.
     * @param kmer El k-mer codificado (uint64_t).
     * @return La frecuencia estimada (CounterType).
     */
    CounterType estimate(uint64_t kmer) const {
        std::vector<CounterType> estimates;
        estimates.reserve(W);

        for (int i = 0; i < W; ++i) {
            uint64_t hash_h = fast_hash(kmer, seeds_h[i]);
            int column_index = hash_h % D;

            // La estimación es el valor del contador en esa posición
            estimates.push_back(matrix[i][column_index]);
        }

        // Devolver la mediana de las estimaciones
        std::sort(estimates.begin(), estimates.end());
        return estimates[W / 2];
    }

    // Métodos para obtener parámetros para el cálculo del Z-Score (en la siguiente etapa)
    int getW() const { return W; }
    int getD() const { return D; }
    
    // NOTA: No se expone la matriz para mantener el encapsulamiento.
    // Los métodos para calcular µk y σk requerirán iterar sobre las celdas o un conjunto
    // de k-mers. Por ahora, solo Estimate() es suficiente.
};

// Ejemplo de uso (ejecutable en un prototipo C++):
/*
int main() {
    // Configuración sugerida: w=5, d=2^20 (1048576)
    CountSketch cs(5, 1048576); 

    // K-mer de ejemplo (simulando un k-mer de 31 bases codificado)
    uint64_t kmer1 = 1234567890ULL;
    uint64_t kmer2 = 1234567891ULL; // Muy similar al kmer1 (colisión de hash probable)
    uint64_t kmer3 = 9876543210ULL;

    cs.update(kmer1);
    cs.update(kmer1);
    cs.update(kmer1);
    cs.update(kmer2); // Se espera una frecuencia real de 1

    std::cout << "Frecuencia estimada para kmer1 (real=3): " << cs.estimate(kmer1) << std::endl;
    std::cout << "Frecuencia estimada para kmer2 (real=1): " << cs.estimate(kmer2) << std::endl;
    std::cout << "Frecuencia estimada para kmer3 (real=0): " << cs.estimate(kmer3) << std::endl;

    return 0;
}
*/