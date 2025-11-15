#include "countsketch.cpp"
#include "lector.cpp"
#include "utils.cpp"
#include <filesystem>

class multi_countsketch {
    private:
        std::vector<CountSketch> multi;
        std::vector<int> K_S;
        std::vector<std::string> dataset_files;
        std::string archivo_actual;
        int N;


    public:
    multi_countsketch(int n, const int k_s[]) : N(n){
        multi.resize(N);

        K_S.assign(k_s, k_s + N);
        
        for (int i = 0; i < N; i++){
            CountSketch temp((pow(2,20)), 5);
            multi.emplace_back(temp);
        }

        try {
            for (const auto& entry : std::filesystem::recursive_directory_iterator("datasets")) {
                if (entry.is_regular_file()) dataset_files.push_back(entry.path().string());
            }
        } catch (const std::exception &e) {
            std::cerr << "Error scanning datasets: " << e.what() << std::endl;
            return;
        }
    }

    std::string sgte_archivo(){
        if (dataset_files.empty()) { std::cerr << "No dataset files found" << std::endl; return 0; }
        archivo_actual = dataset_files.back();
        dataset_files.pop_back();

        lectordatasets lector(archivo_actual);
        std::string texto;
        try {
            texto = lector.leerTexto();
        } catch (const std::exception &e) {
            std::cerr << "Error reading "<<archivo_actual<<": "<<e.what()<<"\n";
            return 0;
        }
        return texto;

    }

    void update(std::string& secuencia){

        for (int i = 0; i < N; ++i) {
            int k = K_S[i];
            // ...
            for (size_t j = 0; j <= secuencia.length() - k; ++j) {
                
                std::string_view kmer_str = std::string_view(secuencia).substr(j, k);
                
                // Llama a la función que calcula el k-mer canónico
                uint64_t encoded_kmer = encode_kmer(kmer_str); 
                
                // Actualizar UNICAMENTE el CountSketch para la longitud 'k' actual
                multi[i].update(encoded_kmer); 
            }
        }
    }


};