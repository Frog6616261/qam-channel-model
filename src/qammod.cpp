#include "qammod.hpp"
#include <memory>


Modulator& Modulator::operator=(Modulator other)
{
    swap(other);
    return *this;
}

void Modulator::swap(Modulator& other)
{
    std::swap(this->m_M, other.m_M);
    std::swap(this->m_avg_power, other.m_avg_power);
    std::swap(this->m_int_symbols, other.m_int_symbols);
    std::swap(this->m_qam_symbols, other.m_qam_symbols);
}
