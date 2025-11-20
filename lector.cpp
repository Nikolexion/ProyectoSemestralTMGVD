#include <string>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <cctype>

class lectordatasets{
    private:
        std::string archivo;
    public:
        lectordatasets(const std::string& nombreArchivo) : archivo(nombreArchivo) {}
        std::string getArchivo() const {
            return archivo;
        }
        
        std::string leerTexto() const {
            std::ifstream file(archivo);
            if (!file.is_open()) {
                throw std::runtime_error("No se pudo abrir el archivo: " + archivo);
            }
            
            // 1. Omitir la línea de encabezado
            std::string header_line;
            std::getline(file, header_line);

            // 2. Leer el resto del archivo y filtrar caracteres
            std::stringstream buffer;
            char c;
            while (file.get(c)) {
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
                    buffer << c;
                }
            }
            return buffer.str();
        }
};