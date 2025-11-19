#include "countsketch.cpp"
#include "lector.cpp"
#include "utils.h"
#include <filesystem>
#include <vector>
#include <string>
#include <string_view>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <mutex>
#include <omp.h> // Librería de OpenMP
#include <memory>


class multi_countsketch {
private:
    std::vector<CountSketch> multi;
    std::vector<int> K_S;
    std::vector<std::string> dataset_files;
    std::string archivo_actual;
    int N;
    int W;
    int D;
    // Mutexes: Uno por cada CountSketch para garantizar la seguridad de hilos
    std::vector<std::unique_ptr<std::mutex>> sketch_mutexes;


public:
    multi_countsketch(int n, const int k_s[], int w, int d) : N(n), W(w), D(d) {
        
        // Inicializar K_S (longitudes de k) correctamente
        K_S.assign(k_s, k_s + N);

        for (int i = 0; i < N; i++){
            // Construir CountSketches con (W, D) = (5, 1048576)
            multi.emplace_back(W, D); 
            
            // Inicializar el mutex correspondiente
            sketch_mutexes.emplace_back(std::make_unique<std::mutex>());
        }

        try {
            for (const auto& entry : std::filesystem::recursive_directory_iterator("datasets")) {
                if (entry.is_regular_file()) dataset_files.push_back(entry.path().string());
            }
        } catch (const std::exception &e) {
            std::cerr << "Error scanning datasets: " << e.what() << std::endl;
        }
    }

    /**
     * @brief Retorna la siguiente secuencia del dataset.
     */
    std::string sgte_archivo(){
        if (dataset_files.empty()) { 
            std::cerr << "No dataset files found" << std::endl; 
            return "";
        }
        archivo_actual = dataset_files.back();
        dataset_files.pop_back();

        lectordatasets lector(archivo_actual);
        std::string texto;
        try {
            // Nota: Aquí el lector debe filtrar el formato FASTA (encabezados y saltos de línea)
            texto = lector.leerTexto(); 
        } catch (const std::exception &e) {
            std::cerr << "Error reading "<<archivo_actual<<": "<<e.what()<<"\n";
            return "";
        }
        return texto;
    }

    /**
     * @brief Procesa la secuencia dada, actualizando todos los CountSketches 
     * (uno por cada k) en paralelo.
     * @param secuencia La cadena de ADN/ARN a procesar.
     */
    void update(std::string& secuencia){
        if (secuencia.empty()) return;

        for (int i = 0; i < N; ++i) {
            int k = K_S[i];
            
            // Verificar longitud mínima para evitar underflow
            if (secuencia.length() < k) continue; 
            
            // Obtener referencias locales al objeto y al mutex
            CountSketch& current_sketch = multi[i];
            std::mutex& current_mutex = *(sketch_mutexes[i]);
            size_t seq_len = secuencia.length();
            
            // ======================== PARALELIZACIÓN ========================
            // Dividimos el trabajo sobre las posiciones (j) entre los hilos
            // Usamos schedule(static) para una distribución de carga eficiente y simple
            #pragma omp parallel for shared(current_sketch, current_mutex, secuencia) schedule(static)
            for (size_t j = 0; j <= seq_len - k; ++j) {
                
                // 1. Generar la vista de k-mer
                std::string_view kmer_str = std::string_view(secuencia).substr(j, k);
                
                // 2. Codificar al k-mer canónico (seguro y local al hilo)
                uint64_t encoded_kmer = encode_kmer(kmer_str); 
                
                // 3. SINCRONIZACIÓN: Protegemos la escritura en el CountSketch
                {
                    // Bloquea el mutex al entrar y lo libera al salir del scope
                    std::lock_guard<std::mutex> lock(current_mutex);
                    
                    // Actualizar el sketch correspondiente
                    current_sketch.update(encoded_kmer); 
                }
            }
        }
    }

    /**
     * @brief Estima la frecuencia de un k-mer dado en el CountSketch correspondiente.
     * @param kmer_str El k-mer en forma de cadena.
     * @param index Índice del CountSketch (0 a N-1).
     * @return La frecuencia estimada del k-mer.
     */
    CounterType estimate(const std::string& kmer_str, int index) {
        if (index < 0 || index >= N) {
            throw std::out_of_range("Index out of range in multi_countsketch::estimate");
        }
        uint64_t encoded_kmer = encode_kmer(kmer_str);
        return multi[index].estimate(encoded_kmer);
    }

    /**
     * @brief Metodo wrapper para ejecutar update en el siguiente archivo del dataset, hasta que se acaben los archivos.
     */
    void procesar_archivos() {
        std::string secuencia = sgte_archivo();
        do {
            update(secuencia);
            secuencia = sgte_archivo();
        } while (!secuencia.empty());
        
    }

    /**
     * @brief Calcula el Score(S) basado en la fórmula de Z-Scores sumados.
     * @param secuencia La secuencia S a evaluar (genoma, lectura, etc.)
     * @param weights (Opcional) Vector de pesos w_k para cada k. Si está vacío, se asume 1.0.
     * @return El puntaje total (double).
     */
    double calculate_score(const std::string& secuencia, const std::vector<double>& weights = {}) {
        double total_score = 0.0;

        // Verificar que tengamos pesos para cada k, si no, usar 1.0
        bool use_custom_weights = (weights.size() == N);

        // Iterar sobre cada configuración de k (k in K)
        for (int i = 0; i < N; ++i) {
            int k = K_S[i];
            double w_k = use_custom_weights ? weights[i] : 1.0;
            
            // Si la secuencia es más corta que k, no podemos sacar k-mers
            if (secuencia.length() < k) continue;

            // 1. Obtener mu_k y sigma_k del sketch actual (multi[i])
            // NOTA: Esto recorre toda la matriz del sketch. Si es muy lento para llamar
            // en cada query, considera calcularlo solo una vez después del entrenamiento.
            std::pair<double, double> stats = multi[i].get_distribution_stats();
            double mu_k = stats.first;
            double sigma_k = stats.second;

            // Evitar división por cero si el sketch está vacío o es uniforme
            if (sigma_k == 0.0) sigma_k = 1.0; 

            double sum_z_scores = 0.0;
            int num_kmers = 0;

            // 2. Sumatoria interna: Recorrer todos los x_k en S
            for (size_t j = 0; j <= secuencia.length() - k; ++j) {
                std::string_view kmer_str = std::string_view(secuencia).substr(j, k);
                uint64_t encoded_kmer = encode_kmer(kmer_str); // Usamos tu función de utils

                // Obtener f_hat (estimación de frecuencia)
                CounterType f_hat = multi[i].estimate(encoded_kmer);

                // Calcular término Z: (f_hat - mu) / sigma
                double z_score = (static_cast<double>(f_hat) - mu_k) / sigma_k;
                
                sum_z_scores += z_score;
                num_kmers++;
            }

            // 3. Aplicar normalización 1 / (|S| - k + 1)
            // El término (|S| - k + 1) es exactamente 'num_kmers'
            double average_z_score = (num_kmers > 0) ? (sum_z_scores / num_kmers) : 0.0;

            // 4. Sumar al score total ponderado
            total_score += w_k * average_z_score;
        }

        return total_score;
    }
};