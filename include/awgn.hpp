#pragma once

#include <vector>
#include <complex>
#include <random>
#include <stdexcept>
#include <cmath>

class GaussianNoiseAdder {
    double m_noise_variance;
    double m_sigma;

    std::mt19937 m_gen;
    std::normal_distribution<double> m_dist;

    
public:
    explicit GaussianNoiseAdder(double noise_variance)
        : m_noise_variance(noise_variance)
        , m_sigma(std::sqrt(noise_variance))
        , m_gen(std::random_device{}())
        , m_dist(0.0, m_sigma)
    {
        if (noise_variance < 0.0) {
            throw std::invalid_argument("Noise variance must be non-negative");
        }
    }

    void add_noise(std::vector<std::complex<double>>& signal)
    {
        for (auto& symb : signal) {
            double noise_re = m_dist(m_gen);
            double noise_im = m_dist(m_gen);

            symb += std::complex<double>(noise_re, noise_im);
        }
    }

    std::vector<std::complex<double>> get_noisy_signal(
        const std::vector<std::complex<double>>& signal
    )
    {
        std::vector<std::complex<double>> noisy_signal = signal;
        add_noise(noisy_signal);
        return noisy_signal;
    }
};