#include <string_view>
#include <string>
#include <algorithm>
#include <cstdint>

char get_base_complement(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return 'N'; // Manejar bases desconocidas
    }
}

std::string reverse_complement(std::string_view kmer_str) {
    std::string rc_kmer;
    rc_kmer.resize(kmer_str.length());
    int k = kmer_str.length();
    
    // Recorrer la cadena al revés
    for (int i = 0; i < k; ++i) {
        // Colocar el complemento en la posición inversa (reverso)
        rc_kmer[k - 1 - i] = get_base_complement(kmer_str[i]);
    }
    return rc_kmer;
}

uint64_t encode_kmer(std::string_view kmer_str) {
    // Definimos una función local para la codificación binaria
    auto binary_encode = [](std::string_view s) -> uint64_t {
        uint64_t encoded = 0;
        for (char base : s) {
            encoded <<= 2;
            switch (base) {
                case 'C': encoded |= 1; break; 
                case 'G': encoded |= 2; break;
                case 'T': encoded |= 3; break;
                // 'A' queda como 00
            }
        }
        return encoded;
    };
    uint64_t kmer_code = binary_encode(kmer_str);
    std::string rc_str = reverse_complement(kmer_str);
    uint64_t rc_code = binary_encode(rc_str);

    return std::min(kmer_code, rc_code);
}