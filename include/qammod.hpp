#include <vector>
#include <complex>

class Modulator{
private:
    double                              m_avg_power;
    uint                                m_M;
    std::vector<uint>                   m_int_symbols;
    std::vector<std::complex<double>>   m_qam_symbols;


public:

    //Ctor
    Modulator() = delete;
    Modulator(const Modulator&);
    Modulator(uint M, double average_power);

    //Dtor
    ~Modulator() = default;

    Modulator& operator=(Modulator); 

private:

    void swap(Modulator&);      
};

