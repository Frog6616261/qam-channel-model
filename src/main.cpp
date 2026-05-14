#include "awgn.hpp"
#include "qammod.hpp"
#include "qamdemod.hpp"


int main(){

    std::cout << "4" << std::endl;
    
    Modulator mod4(4, 1);
    Modulator mod16(16, 1);
    Modulator mod64(64, 1);

    std::vector<uint> a4 = {0, 1, 2, 3};

    std::vector<std::complex<double>> b4;

    b4 = mod4.modulate(a4);

    for (size_t i = 0; i < 4; ++i)
        std::cout <<  b4[i].real() << " " << b4[i].imag() << std::endl;   


    std::cout << "16" << std::endl;
    
    std::vector<uint> a16(0);
    for (uint i = 0; i < 16; ++i) a16.push_back(i);

    std::vector<std::complex<double>> b16 = mod16.modulate(a16);

    for (size_t i = 0; i < 16; ++i){
        std::cout <<  b16[i].real() << " " << b16[i].imag() << std::endl;
    }  


    std::cout << "64" << std::endl;

    

    std::vector<uint> a64;
    for (uint i = 0; i < 64; ++i) a64.push_back(i);

    std::vector<std::complex<double>> b64 = mod64.modulate(a64);

    for (size_t i = 0; i < 64; ++i){
        std::cout << b64[i].real() << " " << b64[i].imag() << std::endl;
    }  

    Demodulator d4(mod4);
    Demodulator d16(mod16);
    Demodulator d64(mod64);

    std::vector<double> c4 = d4.demodulate_soft(b4, 1.0, Demodulator::MAX_LOG_APP);
    std::vector<bool> c16 = d16.demodulate_hard(b16, 1.0, Demodulator::MIN_NORM);
    std::vector<bool> c64 = d64.demodulate_hard(b64, 1.0, Demodulator::MIN_NORM);

    for (size_t i = 0; i < 4; ++i)
        std::cout << c4[i*2] << " " << c4[i*2 + 1] << std::endl;
    
    for (size_t i = 0; i < 16; ++i)
        std::cout << c16[i*4] << " " << c16[i*4 + 1] << " " << c16[i*4 + 2] << " " << c16[i*4 + 3] << std::endl;

    for (size_t i = 0; i < 64; ++i)
        std::cout << c64[i*6] << " " << c64[i*6 + 1] << " " << c64[i*6 + 2] << " " << c64[i*6 + 3] << " " << c64[i*6 + 4] << " " << c64[i*6 + 5] << std::endl;

    return 0;
}