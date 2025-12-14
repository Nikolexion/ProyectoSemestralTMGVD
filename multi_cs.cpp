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
#include <omp.h>
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

public:
    multi_countsketch(int n, const int k_s[], int w, int d) : N(n), W(w), D(d) {
        
        // Inicializar K_S (longitudes de k)
        K_S.assign(k_s, k_s + N);

        // Construir CountSketches y agregarlos al vector
        for (int i = 0; i < N; i++) multi.emplace_back(W, D); 
        
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
            
            // Obtener referencias locales al objeto 
            CountSketch& current_sketch = multi[i];
            size_t seq_len = secuencia.length();
            
            // Paralelizar el procesamiento de k-mers
            #pragma omp parallel for schedule(static)
            for (size_t j = 0; j <= seq_len - k; ++j) {

                std::string_view kmer_str = std::string_view(secuencia).substr(j, k);
                uint64_t encoded_kmer = encode_kmer(kmer_str); 
                
                current_sketch.update(encoded_kmer); 
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
            std::string().swap(secuencia);
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

        // Verificar que tengamos pesos para cada k, si no, usar default 1.0
        bool use_custom_weights = (weights.size() == N);

        // Iterar sobre cada k
        for (int i = 0; i < N; ++i) {
            int k = K_S[i];
            double w_k = use_custom_weights ? weights[i] : 1.0;
            
            // Si la secuencia es más corta que k, no podemos sacar k-mers y continuamos
            if (secuencia.length() < k) continue;

            // 1. Obtener mu_k y sigma_k del sketch actual (multi[i])
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
                uint64_t encoded_kmer = encode_kmer(kmer_str);

                CounterType f_hat = multi[i].estimate(encoded_kmer);
                double z_score = (static_cast<double>(f_hat) - mu_k) / sigma_k;
                
                sum_z_scores += z_score;
                num_kmers++;
            }

            // 3. Aplicar normalización 1 / (|S| - k + 1)
            double average_z_score = (num_kmers > 0) ? (sum_z_scores / num_kmers) : 0.0;

            // 4. Sumar al score total ponderado
            total_score += w_k * average_z_score;
        }

        return total_score;
    }

    /**
     * @brief Guarda toda la estructura en un archivo .bin
     */
    void save_structure(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) {
            throw std::runtime_error("No se pudo abrir el archivo para escribir: " + filename);
        }

        out.write(reinterpret_cast<const char*>(&N), sizeof(N));
        out.write(reinterpret_cast<const char*>(&W), sizeof(W));
        out.write(reinterpret_cast<const char*>(&D), sizeof(D));

        size_t k_size = K_S.size();
        out.write(reinterpret_cast<const char*>(&k_size), sizeof(k_size));
        out.write(reinterpret_cast<const char*>(K_S.data()), k_size * sizeof(int));

        for (const auto& sketch : multi) {
            sketch.save(out);
        }

        out.close();
        std::cout << "Se guardo la estructura en " << filename << std::endl;
    }

    /**
     * @brief Carga la estructura desde un archivo .bin
     */
    void load_structure(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        if (!in.is_open()) {
            throw std::runtime_error("No se pudo abrir el archivo para leer: " + filename);
        }

        int file_N, file_W, file_D;
        in.read(reinterpret_cast<char*>(&file_N), sizeof(file_N));
        in.read(reinterpret_cast<char*>(&file_W), sizeof(file_W));
        in.read(reinterpret_cast<char*>(&file_D), sizeof(file_D));

        if (file_N != N || file_W != W || file_D != D) {
            throw std::runtime_error("Configuracion incompatible entre archivo y codigo.");
        }

        size_t k_size;
        in.read(reinterpret_cast<char*>(&k_size), sizeof(k_size));
        std::vector<int> file_K_S(k_size);
        in.read(reinterpret_cast<char*>(file_K_S.data()), k_size * sizeof(int));
        
        // Verificar que estamos usando los mismos K
        if (file_K_S != K_S) {
            throw std::runtime_error("Los valores de K del archivo no coinciden con la configuración actual.");
        }

        multi.clear();

        for (int i = 0; i < N; ++i) {
            multi.emplace_back(W, D);
            multi.back().load(in); 
        }
        in.close();
        std::cout << "Se cargo la estructura desde " << filename << std::endl;
    }
};