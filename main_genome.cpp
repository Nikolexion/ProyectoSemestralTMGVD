#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <chrono>
#include <fstream>
#include "multi_cs.cpp"

namespace fs = std::filesystem;

// Configuracion por defecto
int D = 1 << 26; 
int W = 5;
std::vector<int> k_values = {15, 21, 31};

// Archivos
const std::string STRUCTURE_FILE = "multi_countsketch_human_genome.bin";
const std::string DATASET_FOLDER = "datasets";
const std::string CSV_OUTPUT_DIR = "plots/csv";
const std::string CSV_FILENAME = "resultados_scores.csv";

// Elimina espacios y llaves {} de un string
std::string clean_string(std::string s) {
    s.erase(std::remove(s.begin(), s.end(), '{'), s.end());
    s.erase(std::remove(s.begin(), s.end(), '}'), s.end());
    s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
    return s;
}

// Convierte la lista limpia en un vector de int
std::vector<int> parse_int_list(const std::string& input) {
    std::vector<int> values;
    std::stringstream ss(clean_string(input));
    std::string segment;
    while (std::getline(ss, segment, ',')) {
        try {
            values.push_back(std::stoi(segment));
        } catch (...) {
            std::cerr << "Error parseando entero: " << segment << std::endl;
        }
    }
    return values;
}

// Convierte la lista limpia en un vector de double
std::vector<double> parse_double_list(const std::string& input) {
    std::vector<double> values;
    std::stringstream ss(clean_string(input));
    std::string segment;
    while (std::getline(ss, segment, ',')) {
        try {
            values.push_back(std::stod(segment));
        } catch (...) {
            std::cerr << "Error parseando double: " << segment << std::endl;
        }
    }
    return values;
}

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

void print_usage(const char* progName) {
    std::cout << "Uso: " << progName << " <modo> [opciones]\n"
              << "Modos:\n"
              << "  train, score, both\n"
              << "Opciones Requeridas:\n"
              << "  -k {k1,k2...}   Lista de k-mers (ej: 15,21,31)\n"
              << "  -d <num>        Dimension D para el sketch (columnas, ej: 67108864)\n"
              << "  -w <num>        Ancho W para el sketch (filas/hashes, ej: 5)\n"
              << "Opciones Opcionales:\n"
              << "  -p {p1,p2...}   Pesos para scoring (ej: 1.0,1.0,1.5). Default: todos 1.0\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    if (mode != "count" && mode != "score" && mode != "both") {
        std::cerr << "Error: Modo desconocido '" << mode << "'\n";
        print_usage(argv[0]);
        return 1;
    }

    // Argumentos
    std::vector<double> pesos = {};
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (i + 1 < argc) {
            if (arg == "-k") {
                k_values = parse_int_list(argv[++i]);
            } else if (arg == "-d") {
                D = std::stoi(argv[++i]);
            } else if (arg == "-w") {
                W = std::stoi(argv[++i]);
            } else if (arg == "-p") {
                pesos = parse_double_list(argv[++i]);
            }
        }
    }

    // validacion
    if (k_values.empty() || D == 0 || W == 0) {
            std::cerr << "Error: Debes especificar -k, -d, -w para inicializar la estructura antes de cargarla." << std::endl;
            return 1;
    }

    multi_countsketch mcs(k_values.size(), k_values.data(), W, D);

    std::vector<std::string> archivos = obtener_archivos(DATASET_FOLDER);

    if (archivos.empty()) {
        std::cerr << "No se encontraron archivos .fa en " << DATASET_FOLDER << std::endl;
        return 1;
    }

    // Conteo
    if (mode == "conteo" || mode == "both") {
        std::cout << "Iniciando conteo" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();

        mcs.procesar_archivos(); 

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Conteo completado en " << elapsed.count() << " segundos." << std::endl;

        // Guardar estructura
        mcs.save_structure(STRUCTURE_FILE);
    }

    if (mode == "score" || mode == "both") {
        std::cout << "Calculando puntajes" << std::endl;

        // Cargar estructura
        if (mode == "score") {
            if (!fs::exists(STRUCTURE_FILE)) {
                std::cerr << "Error: No se encuentra el archivo " << STRUCTURE_FILE 
                          << ". Ejecuta en modo 'train' o 'both' primero." << std::endl;
                return 1;
            }
            mcs.load_structure(STRUCTURE_FILE);
        }

        // Preparar CSV
        fs::create_directories(CSV_OUTPUT_DIR); // Crea carpetas si no existen
        std::string csv_path = CSV_OUTPUT_DIR + "/" + CSV_FILENAME;
        std::ofstream csvFile(csv_path);
        
        if (!csvFile.is_open()) {
            std::cerr << "Error al crear el archivo CSV en " << csv_path << std::endl;
            return 1;
        }

        // Cabecera del CSV
        csvFile << "Archivo,Score" << std::endl;

        auto start = std::chrono::high_resolution_clock::now();

        int total_files = archivos.size();
        int processed_count = 0;

        std::cout << "Procesando " << total_files << " archivos" << std::endl;

        for (const auto& path : archivos) {
            lectordatasets lector(path);
            std::string secuencia = lector.leerTexto();
            double score = mcs.calculate_score(secuencia, pesos);
            
            // Guardar en CSV
            std::string filename = fs::path(path).filename().string();
            csvFile << filename << "," << score << std::endl;

            // Barra de progreso visual
            processed_count++;
            double progress = (double)processed_count / total_files * 100.0;
            std::cout << "\r[" << processed_count << "/" << total_files << "] " 
                      << std::fixed << std::setprecision(1) << progress << "% completado " 
                      << "- Procesando: " << filename << std::string(10, ' ') << std::flush;
        }
        std::cout << std::endl;

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        
        csvFile.close();
        std::cout << "Scoring completado en " << elapsed.count() << " segundos." << std::endl;
        std::cout << "Resultados guardados en: " << csv_path << std::endl;
    }

    return 0;
}