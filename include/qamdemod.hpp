#pragma once
#include "qammod.hpp"


class Demodulator: public Modulator{
private:
    std::vector<std::vector<uint>> m_num_symb_by_bit0_in_pos; // [position bit0 in symb][num qam symbol wit bit0 in position]
    std::vector<std::vector<uint>> m_num_symb_by_bit1_in_pos; // [position bit1 in symb][num qam symbol wit bit0 in position]

public:
    enum MethodsOfDemodulate{
        MIN_NORM,
        SUM_EXP,
        MAX_LOG_APP
    };


public:
    Demodulator() = delete;
    Demodulator(const Demodulator&) = default;
    explicit Demodulator(uint M, double average_power);
    explicit Demodulator(const std::vector<uint>& bit_symbols, const std::vector<std::complex<double>>& qam_symbols);
    explicit Demodulator(const Modulator&);
    Demodulator& operator=(Demodulator);

    ~Demodulator() = default;

    std::vector<bool> demodulate_hard(const std::vector<std::complex<double>>& qam_symbols, double noize_disp, MethodsOfDemodulate);
    std::vector<double> demodulate_soft(const std::vector<std::complex<double>>& qam_symbols, double noize_disp, MethodsOfDemodulate);

private:

    double solve_llr_by_sum_exp(uint bit_num, const std::complex<double>& rx_symb, double noise_disp);
    double solve_llr_by_max_log(uint bit_num, const std::complex<double>& rx_symb, double noise_disp);
    uint find_symb_by_min_norm(const std::complex<double>& rx_symb);

    void init_num_bits_in_symb_arr();

    void swap(Demodulator&);
};

