#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <random>
#include <stdexcept>

#include "qammod.hpp"
#include "qamdemod.hpp"
#include "awgn.hpp"

double convert_dB_to_lin(double dB);

double convert_lin_to_dB(double lin);

double solve_ber(const std::vector<bool>& rx, const std::vector<bool>& tx);

std::vector<bool> generate_random_bits(std::size_t bit_count);

double calculate_noise_sigma(double EbN0_dB, unsigned int bits_per_symbol, double E_avg);

std::size_t count_bit_errors(const std::vector<bool>& rx, const std::vector<bool>& tx);

double simulate_ber_for_one_ebn0(
    Modulator& mod,
    Demodulator& demod,
    double ebn0_db,
    unsigned int bits_per_symbol,
    double E_avg,
    std::size_t symbols_per_block,
    std::size_t min_errors,
    std::size_t max_bits,
    std::size_t& total_errors,
    std::size_t& total_bits
);



int main()
{
    std::vector<double> EbN0_dB;

    for (size_t i = 0; i < 11; ++i)
    {
        EbN0_dB.push_back(static_cast<double>(i*2));
    }

    double E_avg = 1.0;

    Modulator m4(4, E_avg);
    Modulator m16(16, E_avg);
    Modulator m64(64, E_avg);

    Demodulator d4(m4);
    Demodulator d16(m16);
    Demodulator d64(m64);

    const unsigned int bits_per_symbol_4 = 2;
    const unsigned int bits_per_symbol_16 = 4;
    const unsigned int bits_per_symbol_64 = 6;

    const std::size_t symbols_per_block = 100000;

    const std::size_t min_errors = 100;
    const std::size_t max_bits = 10000000; 

    std::ofstream fout("ber_result.txt");

    if (!fout.is_open())
    {
        std::cerr << "Cannot open output file" << std::endl;
        return 1;
    }

    fout << "EbN0_dB "
        << "BER_QAM4 ERR_QAM4 BITS_QAM4 "
        << "BER_QAM16 ERR_QAM16 BITS_QAM16 "
        << "BER_QAM64 ERR_QAM64 BITS_QAM64\n";

    for (double ebn0_db : EbN0_dB)
    {
        size_t errors4 = 0;
        size_t bits4 = 0;

        size_t errors16 = 0;
        size_t bits16 = 0;

        size_t errors64 = 0;
        size_t bits64 = 0;

        const double ber4 = simulate_ber_for_one_ebn0(
            m4,
            d4,
            ebn0_db,
            bits_per_symbol_4,
            E_avg,
            symbols_per_block,
            min_errors,
            max_bits,
            errors4,
            bits4
        );

        const double ber16 = simulate_ber_for_one_ebn0(
            m16,
            d16,
            ebn0_db,
            bits_per_symbol_16,
            E_avg,
            symbols_per_block,
            min_errors,
            max_bits,
            errors16,
            bits16
        );

        const double ber64 = simulate_ber_for_one_ebn0(
            m64,
            d64,
            ebn0_db,
            bits_per_symbol_64,
            E_avg,
            symbols_per_block,
            min_errors,
            max_bits,
            errors64,
            bits64
        );

        std::cout
            << "Eb/N0 = " << ebn0_db << " dB"
            << " | QAM-4 BER = " << ber4
            << " errors = " << errors4
            << " bits = " << bits4
            << " | QAM-16 BER = " << ber16
            << " errors = " << errors16
            << " bits = " << bits16
            << " | QAM-64 BER = " << ber64
            << " errors = " << errors64
            << " bits = " << bits64
            << std::endl;

        fout
            << ebn0_db << " "
            << ber4 << " " << errors4 << " " << bits4 << " "
            << ber16 << " " << errors16 << " " << bits16 << " "
            << ber64 << " " << errors64 << " " << bits64 << "\n";
    }

    fout.close();

    std::cout << "Results saved to ber_result.txt" << std::endl;

    return 0;
}




double convert_dB_to_lin(double dB)
{
    return std::pow(10.0, dB / 10.0);
}

double convert_lin_to_dB(double lin)
{
    if (lin <= 0.0)
    {
        throw std::invalid_argument("lin must be positive");
    }

    return 10.0 * std::log10(lin);
}

double solve_ber(const std::vector<bool>& rx, const std::vector<bool>& tx)
{
    if (rx.size() != tx.size()) throw std::invalid_argument("rx and tx must have the same size");


    if (tx.empty()) throw std::invalid_argument("rx and tx must not be empty");

    std::size_t errors = 0;

    for (std::size_t i = 0; i < tx.size(); ++i){
        if (rx[i] != tx[i]){
            ++errors;
        }
    }

    return static_cast<double>(errors) / static_cast<double>(tx.size());
}

std::vector<bool> generate_random_bits(std::size_t bit_count)
{
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_int_distribution<int> dist(0, 1);

    std::vector<bool> bits;
    bits.reserve(bit_count);

    for (std::size_t i = 0; i < bit_count; ++i){
        bits.push_back(static_cast<bool>(dist(gen)));
    }

    return bits;
}

double calculate_noise_sigma(double EbN0_dB, unsigned int bits_per_symbol, double E_avg)
{
    const double EbN0_lin = convert_dB_to_lin(EbN0_dB);

    const double Eb = E_avg / static_cast<double>(bits_per_symbol);

    const double N0 = Eb / EbN0_lin;

    const double sigma2 = N0 / 2.0;
    const double sigma = std::sqrt(sigma2);

    return sigma;
}

double calculate_noise_variance(double EbN0_dB, unsigned int bits_per_symbol, double E_avg)
{
    const double EbN0_lin = convert_dB_to_lin(EbN0_dB);

    const double Eb = E_avg / static_cast<double>(bits_per_symbol);

    const double N0 = Eb / EbN0_lin;

    const double noise_variance = N0 / 2.0;

    return noise_variance;
}

std::size_t count_bit_errors(const std::vector<bool>& rx, const std::vector<bool>& tx)
{
    if (rx.size() != tx.size()) throw std::invalid_argument("rx and tx must have the same size");

    size_t errors = 0;

    for (size_t i = 0; i < tx.size(); ++i)
    {
        if (rx[i] != tx[i])
        {
            ++errors;
        }
    }

    return errors;
}

double simulate_ber_for_one_ebn0(
    Modulator& mod,
    Demodulator& demod,
    double ebn0_db,
    unsigned int bits_per_symbol,
    double E_avg,
    size_t symbols_per_block,
    size_t min_errors,
    size_t max_bits,
    size_t& total_errors,
    size_t& total_bits
)
{
    total_errors = 0;
    total_bits = 0;

    const double noise_var = calculate_noise_variance(ebn0_db, bits_per_symbol, E_avg);



    while (total_errors < min_errors && total_bits < max_bits)
    {
        const size_t bits_in_block = symbols_per_block * bits_per_symbol;

        std::vector<bool> tx = generate_random_bits(bits_in_block);

        std::vector<std::complex<double>> symbols = mod.modulate(tx);

        GaussianNoiseAdder awgn(noise_var);
        awgn.add_noise(symbols);

        std::vector<bool> rx = demod.demodulate_hard(
            symbols,
            noise_var,
            Demodulator::MIN_NORM
        );

        total_errors += count_bit_errors(rx, tx);
        total_bits += tx.size();
    }

    if (total_bits == 0)
    {
        return 0.0;
    }

    return static_cast<double>(total_errors) / static_cast<double>(total_bits);
}