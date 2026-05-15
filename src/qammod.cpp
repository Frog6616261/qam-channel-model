#include "qammod.hpp"
#include <memory>
#include <cmath>
#include <numeric>
#include <stdexcept>  
#include <algorithm>

double meanEnergy(const std::vector<std::complex<double>>&);
bool isPowerOfTwo(uint);

uint log2_PowerOfTwo(uint value);

template<typename T>
bool allDifferent(std::vector<T>);


Modulator::Modulator(uint M, double average_power)
    : m_M(M)
    , m_avg_power(average_power)
    , m_num_qam_by_bits_symbols(M)
    , m_qam_symbols(M)
{
    init_constellation_map(m_M, m_avg_power);
}

Modulator::Modulator(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& complex_symbols)
    : m_M(bit_symbols.size())
    , m_avg_power(0.0)
    , m_num_qam_by_bits_symbols(bit_symbols)
    , m_qam_symbols(complex_symbols)
{
    init_constellation_map(bit_symbols, complex_symbols);
    m_avg_power = meanEnergy(complex_symbols);
}


Modulator& Modulator::operator=(Modulator other)
{
    swap(other);
    return *this;
}

void Modulator::set_M(uint M)
{
    m_M = M;
    m_num_qam_by_bits_symbols.assign(M, 0);
    m_qam_symbols.assign(M, {0.0, 0.0});
}

void Modulator::set_avg_power(double avg_power)
{
    init_constellation_map(m_M, avg_power);
}

uint Modulator::get_M()
{
    return m_M;
}

double Modulator::get_avg_power()
{
    return m_avg_power;
}

std::vector<uint> Modulator::get_int_symbols()
{
    return m_num_qam_by_bits_symbols;
}

std::vector<std::complex<double>> Modulator::get_qam_symbols()
{
    return m_qam_symbols;
}


std::vector<std::complex<double>> Modulator::modulate(std::vector<uint> int_symbols)
{
    uint sz = int_symbols.size();
    std::vector<std::complex<double>> result(sz, {0.0, 0.0});

    for (uint num_int_symb = 0; num_int_symb < sz; ++num_int_symb){
        result[num_int_symb] = m_qam_symbols[m_num_qam_by_bits_symbols[int_symbols[num_int_symb]]];
    }

    return result;
}

std::vector<std::complex<double>> Modulator::modulate(std::vector<bool> bits)
{
    const uint bits_per_symbol = log2_PowerOfTwo(m_qam_symbols.size());

    if (bits.size() % bits_per_symbol != 0){
        throw std::invalid_argument("Number of bits must be divisible by bits_per_symbol");
    }

    const uint symbols_count = static_cast<uint>(bits.size() / bits_per_symbol);

    std::vector<uint> int_symbols(symbols_count, 0);

    for (uint i = 0; i < symbols_count; ++i)
    {
        uint value = 0;

        for (uint bit_num = 0; bit_num < bits_per_symbol; ++bit_num)
        {
            value <<= 1;

            if (bits[i * bits_per_symbol + bit_num])
            {
                value |= 1u;
            }
        }

        int_symbols[i] = value;
    }

    return modulate(int_symbols);
}


void Modulator::init_constellation_map(uint M, double average_power)
{
    if (!isPowerOfTwo(M))           throw std::invalid_argument("M must be power of 2");
    if (average_power <= 0)         throw std::invalid_argument("Average Power must be positive");

    switch (M)
    {
    case 4:
        fill_symbols_by_qam4();
        break;
    
    case 16:
        fill_symbols_by_qam16();
        break;

    case 64:
        fill_symbols_by_qam64();
        break;

    default:
        throw std::invalid_argument("M must be 4, 16 or 64");
        break;
    }
}

void Modulator::init_constellation_map(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& complex_symbols)
{
    if (bit_symbols.size() != complex_symbols.size())                   throw std::logic_error("Vectors are not equal");
    if (!isPowerOfTwo(bit_symbols.size()))                              throw std::invalid_argument("M must be power of 2");
    if (meanEnergy(complex_symbols) <= 0)                               throw std::invalid_argument("Average Power must be positive");
    if (!allDifferent(bit_symbols) || !allDifferent(complex_symbols))     throw std::logic_error("Vector's has not different elements");

    m_num_qam_by_bits_symbols = bit_symbols;
    m_qam_symbols = complex_symbols;
}

void Modulator::fill_symbols_by_qam4()
{
    for (uint int_numb = 0; int_numb < 4; ++int_numb){
        double bit0_val = static_cast<double>(int_numb & 1);
        double bit1_val = static_cast<double>((int_numb >> 1) & 1);

        m_num_qam_by_bits_symbols[int_numb] = int_numb;
        m_qam_symbols[int_numb] = (std::sqrt(m_avg_power)/std::sqrt(2.0)) 
        * std::complex<double>((1.0 - 2.0*bit0_val), (1.0 - 2.0*bit1_val));
    }
}

void Modulator::fill_symbols_by_qam16()
{
    for (uint int_numb = 0; int_numb < 16; ++int_numb){
        double bit0_val = static_cast<double>(int_numb & 1);
        double bit1_val = static_cast<double>((int_numb >> 1) & 1);
        double bit2_val = static_cast<double>((int_numb >> 2) & 1);
        double bit3_val = static_cast<double>((int_numb >> 3) & 1);

        m_num_qam_by_bits_symbols[int_numb] = int_numb;
        m_qam_symbols[int_numb] = (std::sqrt(m_avg_power)/std::sqrt(10.0)) 
        * std::complex<double>((1.0 - 2.0*bit0_val)*(2.0 - (1.0 - 2.0*bit2_val))
        , (1.0 - 2.0*bit1_val)*(2.0 - (1.0 - 2.0*bit3_val)));
    }
}

void Modulator::fill_symbols_by_qam64()
{
    for (uint int_numb = 0; int_numb < 64; ++int_numb){
        double bit0_val = static_cast<double>(int_numb & 1);
        double bit1_val = static_cast<double>((int_numb >> 1) & 1);
        double bit2_val = static_cast<double>((int_numb >> 2) & 1);
        double bit3_val = static_cast<double>((int_numb >> 3) & 1);
        double bit4_val = static_cast<double>((int_numb >> 4) & 1);
        double bit5_val = static_cast<double>((int_numb >> 5) & 1);

        m_num_qam_by_bits_symbols[int_numb] = int_numb;
        m_qam_symbols[int_numb] = (std::sqrt(m_avg_power)/std::sqrt(42.0)) 
        * std::complex<double>((1.0 - 2.0*bit0_val)*(4.0 - (1.0 - 2.0*bit2_val)*(2.0 - (1.0 - 2.0*bit4_val)))
        , (1.0 - 2.0*bit1_val)*(4.0 - (1.0 - 2.0*bit3_val)*(2.0 - (1.0 - 2.0*bit5_val))));
    }
}

void Modulator::swap(Modulator& other)
{
    std::swap(m_M, other.m_M);
    std::swap(m_avg_power, other.m_avg_power);
    std::swap(m_num_qam_by_bits_symbols, other.m_num_qam_by_bits_symbols);
    std::swap(m_qam_symbols, other.m_qam_symbols);
}




double meanEnergy(const std::vector<std::complex<double>>& symbols)
{
    if (symbols.empty()) {
        throw std::invalid_argument("Vector is empty");
    }

    return std::accumulate(symbols.begin(), symbols.end(), 0.0,
        [](double sum, const std::complex<double>& x) {return sum + std::norm(x);})
        / symbols.size();
}

bool isPowerOfTwo(uint n)
{
    return n > 0 && (n & (n - 1)) == 0;
}

uint log2_PowerOfTwo(uint value)
{
    if (value == 0 || (value & (value - 1)) != 0) {
        throw std::invalid_argument("Value must be a positive power of two");
    }

    uint result = 0;

    while (value > 1) {
        value >>= 1;
        ++result;
    }

    return result;
}

template<typename T>
bool allDifferent(std::vector<T> v) {
    for (size_t i = 0; i < v.size(); ++i) {
        for (size_t j = i + 1; j < v.size(); ++j) {
            if (v[i] == v[j]) {
                return false;
            }
        }
    }

    return true;
}