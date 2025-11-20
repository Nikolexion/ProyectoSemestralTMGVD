#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include "multi_cs.cpp"

namespace fs = std::filesystem;

// --- CONFIGURACIÓN PARA GENOMA HUMANO ---
// 2^26 = ~67 millones de columnas. 
// Memoria aprox: 67M * 8 bytes * 5 filas = ~2.6 GB de RAM. 
// Bajar a 25 (1.3 GB) o 24 (600 MB) según la RAM disponible.
const int HUMAN_D = 1 << 26; 
const int HUMAN_W = 5;

// Función auxiliar para obtener todos los archivos .fa
std::vector<std::string> obtener_archivos(const std::string& ruta) {
    std::vector<std::string> archivos;
    try {
        for (const auto& entry : fs::directory_iterator(ruta)) {
            if (entry.path().extension() == ".fa" || entry.path().extension() == ".fasta") {
                archivos.push_back(entry.path().string());
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error leyendo directorio: " << e.what() << std::endl;
    }
    return archivos;
}

int main() {
    // 1. Configuración de K-mers (ej. k=21 para especificidad en genomas grandes)
    // Probamos con las combinaciones {15, 21, 31}
    const int N = 3;
    int k_values[N] = {15, 21, 31};

    std::cout << "=== PROCESAMIENTO GENOMA HUMANO GRCh38 ===" << std::endl;
    std::cout << "Configuracion: W=" << HUMAN_W << ", D=" << HUMAN_D << std::endl;
    
    // Instanciar el Sketch
    multi_countsketch mcs(N, k_values, HUMAN_W, HUMAN_D);

    std::string carpeta_datasets = "datasets";
    std::vector<std::string> archivos = obtener_archivos(carpeta_datasets);

    if (archivos.empty()) {
        std::cerr << "No se encontraron archivos .fa en " << carpeta_datasets << std::endl;
        return 1;
    }

    // ==========================================
    // FASE 1: CONTEO APROXIMADO DE K-MERS
    // ==========================================
    std::cout << "\n--- FASE 1: CONTEO APROXIMADO DE K-MERS ---" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    mcs.procesar_archivos();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Conteo completado en " << elapsed.count() << " segundos." << std::endl;


    // ==========================================
    // FASE 2: PUNTUACIÓN (Scoring)
    // ==========================================
    
    std::cout << "\n--- FASE 2: CALCULO DE SCORES ---" << std::endl;
    std::cout << "Archivo | Score (Indice de Anomalia/Rareza)" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    // Pesos opcionales (dando mas peso a k-mers largos)
    std::vector<double> pesos = {1.0, 1.5, 2.0}; 

    for (const auto& path : archivos) {
        lectordatasets lector(path);
        std::string secuencia = lector.leerTexto();
        double score = mcs.calculate_score(secuencia, pesos);
        std::cout << fs::path(path).filename().string() << " | " << score << std::endl;
    }
    std::cout << "\nProceso finalizado." << std::endl;
    return 0;
}