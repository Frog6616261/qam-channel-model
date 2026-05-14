#pragma once

#include <vector>
#include <complex>

class Modulator{
protected:
    uint                                m_M;
    double                              m_avg_power;
    std::vector<uint>                   m_int_symbols;
    std::vector<std::complex<double>>   m_qam_symbols;


public:

    //Ctor
    Modulator() = delete;
    Modulator(const Modulator&) = default;
    explicit Modulator(uint M, double average_power);
    explicit Modulator(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& complex_symbols);

    //Dtor
    ~Modulator() = default;

    Modulator& operator=(Modulator); 

    void set_M(uint M);
    void set_avg_power(double avg_power);

    uint get_M();
    double get_avg_power();
    std::vector<uint> get_int_symbols();
    std::vector<std::complex<double>> get_qam_symbols();

    void set_castom_constellation_map(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& complex_symbols);
    std::vector<std::complex<double>> modulate(std::vector<uint>);

private:

    void init_constellation_map(uint M, double average_power);
    void init_constellation_map(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& complex_symbols);
    
    void fill_symbols_by_qam4();
    void fill_symbols_by_qam16();
    void fill_symbols_by_qam64();

    void swap(Modulator&);      
};

