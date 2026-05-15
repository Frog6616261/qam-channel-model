import sys
import numpy as np
import matplotlib.pyplot as plt


def replace_zeros_for_log_plot(ber):
    """
    На semilogy нельзя нормально отобразить BER = 0.
    Заменяем нули на NaN, чтобы matplotlib их пропускал.
    """
    ber = np.array(ber, dtype=float)
    ber[ber <= 0.0] = np.nan
    return ber


def plot_ber(filename):
    data = np.loadtxt(filename, skiprows=1)

    ebn0_db = data[:, 0]

    ber_qam4 = data[:, 1]
    err_qam4 = data[:, 2]
    bits_qam4 = data[:, 3]

    ber_qam16 = data[:, 4]
    err_qam16 = data[:, 5]
    bits_qam16 = data[:, 6]

    ber_qam64 = data[:, 7]
    err_qam64 = data[:, 8]
    bits_qam64 = data[:, 9]

    plt.figure(figsize=(10, 6))

    plt.semilogy(
        ebn0_db,
        replace_zeros_for_log_plot(ber_qam4),
        marker="o",
        linestyle="-",
        label="QAM4"
    )

    plt.semilogy(
        ebn0_db,
        replace_zeros_for_log_plot(ber_qam16),
        marker="s",
        linestyle="-",
        label="QAM16"
    )

    plt.semilogy(
        ebn0_db,
        replace_zeros_for_log_plot(ber_qam64),
        marker="^",
        linestyle="-",
        label="QAM64"
    )

    plt.grid(True, which="both")
    plt.xlabel("Eb/N0, dB")
    plt.ylabel("BER")
    plt.title("BER curves for QAM modulations")
    plt.legend()
    plt.tight_layout()

    plt.savefig("ber_plot.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else "../ber_result.txt"
    plot_ber(filename)