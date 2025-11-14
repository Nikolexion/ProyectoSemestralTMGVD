#include <string>
#include <fstream>
#include <stdexcept>

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
            std::string contenido((std::istreambuf_iterator<char>(file)),
                                 std::istreambuf_iterator<char>());
            file.close();
            return contenido;
        }
};