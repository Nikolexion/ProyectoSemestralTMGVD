#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <random>

// Definimos el tipo de contador
using CounterType = uint32_t;

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
            //matrix[i][column_index]++;
            matrix[i][column_index] += sign;
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
            estimates.push_back(matrix[i][column_index]);
        }

        // Devolver la mediana de las estimaciones
        std::sort(estimates.begin(), estimates.end());
        return estimates[W / 2];
    }

    // Método para obtener el parámetro w 
    int getW() const { return W; }
    // Método para obtener el parámetro d
    int getD() const { return D; }

    /**
     * @brief Calcula la media y la desviación estándar de los contadores en la matriz.
     * Esto sirve para normalizar los puntajes (Z-Score).
     * @return Un par {media, desviacion_estandar}
     */
    std::pair<double, double> get_distribution_stats() const {
        double sum = 0.0;
        double sum_sq = 0.0;
        long long total_elements = (long long)W * D;
        for (const auto& row : matrix) {
            for (CounterType val : row) {
                double v = static_cast<double>(val);
                sum += v;
                sum_sq += (v * v);
            }
        }

        // Calcular media
        double mean = sum / total_elements;

        // Calcular varianza y desviación estándar
        double variance = (sum_sq / total_elements) - (mean * mean);
        
        // Evitar raíces negativas por errores de punto flotante muy pequeños
        if (variance < 0) variance = 0; 

        double std_dev = std::sqrt(variance);

        return {mean, std_dev};
    }
};