#include "qamdemod.hpp"
#include <memory>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <algorithm>


uint log2PowerOfTwo(uint);


Demodulator::Demodulator(uint M, double average_power)
    : Demodulator::Modulator(M, average_power)
    , m_num_symb_by_bit0_in_pos(log2PowerOfTwo(m_M), std::vector<uint>(m_M >> 1))
    , m_num_symb_by_bit1_in_pos(log2PowerOfTwo(m_M), std::vector<uint>(m_M >> 1))
{
    init_num_bits_in_symb_arr();
}

Demodulator::Demodulator(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& qam_symbols)
    : Demodulator::Modulator(bit_symbols, qam_symbols)
    , m_num_symb_by_bit0_in_pos(log2PowerOfTwo(m_M), std::vector<uint>(m_M >> 1))
    , m_num_symb_by_bit1_in_pos(log2PowerOfTwo(m_M), std::vector<uint>(m_M >> 1))
{
    init_num_bits_in_symb_arr();
}

Demodulator::Demodulator(const Modulator& other)
    : Demodulator::Modulator(other)
    , m_num_symb_by_bit0_in_pos(log2PowerOfTwo(m_M), std::vector<uint>(m_M >> 1))
    , m_num_symb_by_bit1_in_pos(log2PowerOfTwo(m_M), std::vector<uint>(m_M >> 1))
{
    init_num_bits_in_symb_arr();
}

Demodulator& Demodulator::operator=(Demodulator other)
{
    swap(other);
    return *this;
}



std::vector<bool> Demodulator::demodulate_hard(const std::vector<std::complex<double>>& qam_symbols
    , double noize_disp, MethodsOfDemodulate demod_method)
{
    if (noize_disp <= 0.0) throw std::invalid_argument("Noise variance must be positive");


    uint sz = qam_symbols.size();
    uint numb_bits_in_symb = log2PowerOfTwo(m_M);
    std::vector<bool> result(sz*numb_bits_in_symb, 0);

    switch (demod_method)
    {
    case MIN_NORM:
        for (uint num_symb = 0; num_symb < sz; ++num_symb){
            uint cur_demod_symb = find_symb_by_min_norm(qam_symbols[num_symb]);

            for (uint num_bit_in_symb = 0; num_bit_in_symb < numb_bits_in_symb; ++num_bit_in_symb){
                uint num_bit = numb_bits_in_symb - num_bit_in_symb - 1;
                result[num_symb*numb_bits_in_symb + num_bit_in_symb] = (cur_demod_symb >> num_bit) & 1;
            }
        }
        break;

    case SUM_EXP:
        for (uint num_symb = 0; num_symb < sz; ++num_symb){
            for (uint num_bit_in_symb = 0; num_bit_in_symb < numb_bits_in_symb; ++num_bit_in_symb){

                result[num_symb*numb_bits_in_symb + num_bit_in_symb] = 
                solve_llr_by_sum_exp(num_bit_in_symb, qam_symbols[num_symb], noize_disp) > 0;
            }
        }
        break;

    case MAX_LOG_APP:
        for (uint num_symb = 0; num_symb < sz; ++num_symb){
            for (uint num_bit_in_symb = 0; num_bit_in_symb < numb_bits_in_symb; ++num_bit_in_symb){

                result[num_symb*numb_bits_in_symb + num_bit_in_symb] = 
                solve_llr_by_max_log(num_bit_in_symb, qam_symbols[num_symb], noize_disp) > 0;
            }
        }
        break;

    default:
        break;
    }

    return result;
}

std::vector<double> Demodulator::demodulate_soft(const std::vector<std::complex<double>>& qam_symbols
    , double noize_disp, MethodsOfDemodulate demod_method)
{
    if (noize_disp <= 0.0) throw std::invalid_argument("Noise variance must be positive");
    
    uint sz = qam_symbols.size();
    uint numb_bits_in_symb = log2PowerOfTwo(m_M);
    std::vector<double> result(sz*numb_bits_in_symb, 0);

    switch (demod_method)
    {
    case SUM_EXP:
        for (uint num_symb = 0; num_symb < sz; ++num_symb){
            for (uint num_bit_in_symb = 0; num_bit_in_symb < numb_bits_in_symb; ++num_bit_in_symb){

                result[num_symb*numb_bits_in_symb + num_bit_in_symb] = 
                solve_llr_by_sum_exp(num_bit_in_symb, qam_symbols[num_symb], noize_disp);
            }
        }
        break;

    case MAX_LOG_APP:
        for (uint num_symb = 0; num_symb < sz; ++num_symb){
            for (uint num_bit_in_symb = 0; num_bit_in_symb < numb_bits_in_symb; ++num_bit_in_symb){

                result[num_symb*numb_bits_in_symb + num_bit_in_symb] = 
                solve_llr_by_max_log(num_bit_in_symb, qam_symbols[num_symb], noize_disp);
            }
        }
        break;

    default:
        break;
    }

    return result;
}


double Demodulator::solve_llr_by_sum_exp(uint num_bit_in_symb
    , const std::complex<double>& rx_symb, double noise_disp)
{
    if (num_bit_in_symb >= log2PowerOfTwo(m_M)) {
        throw std::out_of_range("Bit number is out of range");
    }

    double exp_sum_1 = 0.0;
    double exp_sum_0 = 0.0;
    uint numb_symbs_for_bit = m_M >> 1;

    for (uint num_symb_in_arr = 0; num_symb_in_arr < numb_symbs_for_bit; ++num_symb_in_arr){
        const double norm_bit1 = std::norm(
            m_qam_symbols[m_num_symb_by_bit1_in_pos[num_bit_in_symb][num_symb_in_arr]] - rx_symb);
        const double norm_bit0 = std::norm(
            m_qam_symbols[m_num_symb_by_bit0_in_pos[num_bit_in_symb][num_symb_in_arr]] - rx_symb);

        exp_sum_1 += std::exp(-norm_bit1 / (2.0 * noise_disp));
        exp_sum_0 += std::exp(-norm_bit0 / (2.0 * noise_disp));
    }

    return std::log(exp_sum_1) - std::log(exp_sum_0);
}


double Demodulator::solve_llr_by_max_log(uint num_bit_in_symb
    , const std::complex<double>& rx_symb, double noise_disp)
{
    double min_norm_0 = std::numeric_limits<double>::max();
    double min_norm_1 = std::numeric_limits<double>::max();

    uint numb_symbs_for_bit = m_M >> 1;

    double cur_norm_bit0 = 0;
    double cur_norm_bit1 = 0;
    for (uint num_symb_in_arr = 0; num_symb_in_arr < numb_symbs_for_bit; ++num_symb_in_arr){

        cur_norm_bit0 = std::norm(
            m_qam_symbols[m_num_symb_by_bit0_in_pos[num_bit_in_symb][num_symb_in_arr]] - rx_symb);
        cur_norm_bit1 = std::norm(
            m_qam_symbols[m_num_symb_by_bit1_in_pos[num_bit_in_symb][num_symb_in_arr]] - rx_symb);

        if (cur_norm_bit0 < min_norm_0) min_norm_0 = cur_norm_bit0;
        if (cur_norm_bit1 < min_norm_1) min_norm_1 = cur_norm_bit1;
    }

    return (min_norm_0 - min_norm_1) / (2.0 * noise_disp);
}


uint Demodulator::find_symb_by_min_norm(const std::complex<double>& rx_symb)
{
    std::vector<uint>::iterator it = std::min_element(m_num_qam_by_bits_symbols.begin(), m_num_qam_by_bits_symbols.end(),
    [this, &rx_symb](uint a, uint b) {
        return  std::norm(m_qam_symbols[a] - rx_symb) 
                < std::norm(m_qam_symbols[b] - rx_symb);
    });

    std::vector<uint>::iterator bit_symb_it = std::find(m_num_qam_by_bits_symbols.begin(), m_num_qam_by_bits_symbols.end(), *it);

    return *bit_symb_it;
}


void Demodulator::init_num_bits_in_symb_arr()
{
    uint numb_bits_in_symb = log2PowerOfTwo(m_M);

    for (uint num_bit_in_symb = 0; num_bit_in_symb < numb_bits_in_symb; ++num_bit_in_symb){
        uint cur_num_inbit0_arr = 0;
        uint cur_num_inbit1_arr = 0;

        for (uint num_symb = 0; num_symb < m_M; ++num_symb){

            uint cur_int_symb = m_num_qam_by_bits_symbols[num_symb];
            uint num_bit = numb_bits_in_symb - num_bit_in_symb - 1;
            bool cur_bit_val = (cur_int_symb >> num_bit) & 1;

            if (cur_bit_val){
                m_num_symb_by_bit1_in_pos[num_bit_in_symb][cur_num_inbit1_arr] = num_symb;
                cur_num_inbit1_arr++;
            } else {
                m_num_symb_by_bit0_in_pos[num_bit_in_symb][cur_num_inbit0_arr] = num_symb;
                cur_num_inbit0_arr++;
            }
        }
    }  
}


void Demodulator::swap(Demodulator& other)
{
    std::swap(m_M, other.m_M);
    std::swap(m_avg_power, other.m_avg_power);
    std::swap(m_num_qam_by_bits_symbols, other.m_num_qam_by_bits_symbols);
    std::swap(m_qam_symbols, other.m_qam_symbols);
    std::swap(m_num_symb_by_bit0_in_pos, other.m_num_symb_by_bit0_in_pos);
    std::swap(m_num_symb_by_bit1_in_pos, other.m_num_symb_by_bit1_in_pos);
}


uint log2PowerOfTwo(uint value)
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