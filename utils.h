#ifndef UTILS_CPP
#define UTILS_CPP
#include <string_view>
#include <string>
#include <algorithm>
#include <cstdint>

inline uint64_t base_to_int(char base) {
    switch (base) {
        case 'C': return 1; 
        case 'G': return 2; 
        case 'T': return 3;
        default: return 0; // 'A' y otros
    }
}

uint64_t encode_kmer(std::string_view kmer_str) {
    uint64_t kmer_code = 0;
    uint64_t rc_code = 0;

    // Recorremos el string UNA sola vez para calcular ambos valores
    for (char base : kmer_str) {
        // 1. Codificación Normal (Forward)
        // Desplazamos a la izquierda y agregamos la nueva base
        kmer_code = (kmer_code << 2) | base_to_int(base);
    }

    // 2. Codificación Reverso Complementaria (Reverse Complement)
    // Para el RC, leemos el kmer "al revés" lógicamente.
    // Matemáticamente: el complemento de A(0) es T(3), C(1) es G(2).
    // Nota que 0^3=3, 1^3=2. ¡Es un XOR con 3!
    
    // Iteramos en reverso sobre el string original
    for (auto it = kmer_str.rbegin(); it != kmer_str.rend(); ++it) {
        uint64_t val = base_to_int(*it);
        uint64_t complement_val = val ^ 3; // Truco de bits: XOR 3 da el complemento (0<->3, 1<->2)
        
        rc_code = (rc_code << 2) | complement_val;
    }

    // Retornamos el menor de los dos (Canónico)
    return std::min(kmer_code, rc_code);
}

#endif