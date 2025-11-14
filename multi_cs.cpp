#include "countsketch.cpp"
#include <filesystem>
#include "lector.cpp"

class multi_countsketch {
    private:
        std::vector<CountSketch> multi;
        std::vector<int> K_S;
        std::vector<std::string> dataset_files;
        std::string archivo_actual;


    public:
    multi_countsketch(int n, int k_s[]){
        multi.resize(n);
        K_S.resize(n);

        for (int i = 0; i < n; i++){
            K_S.push_back(k_s[i]);
        }
        
        for (int i = 0; i < n; i++){
            CountSketch temp((pow(2,15)), 3);
            multi.push_back(temp);
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


};