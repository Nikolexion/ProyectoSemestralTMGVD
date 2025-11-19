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

// W: Profundidad (Filas). 5 es un buen valor para precisión.
const int SKETCH_W = 5;       
// D: Ancho (Columnas). 2^20 (1,048,576) es un buen valor para baja colisión.
const int SKETCH_D = 1048576; 

class multi_countsketch {
private:
    std::vector<CountSketch> multi;
    std::vector<int> K_S;
    std::vector<std::string> dataset_files;
    std::string archivo_actual;
    int N;
    // Mutexes: Uno por cada CountSketch para garantizar la seguridad de hilos
    std::vector<std::mutex> sketch_mutexes;


public:
    multi_countsketch(int n, const int k_s[]) : N(n){
        
        // Inicializar K_S (longitudes de k) correctamente
        K_S.assign(k_s, k_s + N);

        for (int i = 0; i < N; i++){
            // Construir CountSketches con (W, D) = (5, 1048576)
            multi.emplace_back(SKETCH_W, SKETCH_D); 
            
            // Inicializar el mutex correspondiente
            sketch_mutexes.emplace_back(); 
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
            std::mutex& current_mutex = sketch_mutexes[i];
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
};