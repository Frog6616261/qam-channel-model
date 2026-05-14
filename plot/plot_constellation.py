import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt


def read_constellations(filename):
    """
    Читает файл формата:

    4
    x1 y1
    x2 y2
    ...

    16
    x1 y1
    ...

    Возвращает список:
    [
        (4, [(x, y), ...]),
        (16, [(x, y), ...]),
        ...
    ]
    """

    constellations = []

    with open(filename, "r", encoding="utf-8") as file:
        lines = [line.strip() for line in file if line.strip()]

    i = 0

    while i < len(lines):
        m = int(lines[i])
        i += 1

        points = []

        for _ in range(m):
            x_str, y_str = lines[i].split()
            x = float(x_str)
            y = float(y_str)

            points.append((x, y))
            i += 1

        constellations.append((m, points))

    return constellations


def bits_for_index(index, m):
    """
    Возвращает битовое представление индекса.

    Например:
    для M = 4 ширина 2 бита:
        0 -> 00
        1 -> 01
        2 -> 10
        3 -> 11

    для M = 16 ширина 4 бита.
    для M = 64 ширина 6 бит.
    """

    bit_count = int(math.log2(m))
    return format(index, f"0{bit_count}b")


def plot_constellation(m, points, output_dir=None, show=True):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]

    plt.figure(figsize=(7, 7))

    plt.scatter(xs, ys, s=60)

    for index, (x, y) in enumerate(points):
        bits = bits_for_index(index, m)

        plt.annotate(
            bits,
            xy=(x, y),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9
        )

    plt.axhline(0, linewidth=1)
    plt.axvline(0, linewidth=1)

    plt.grid(True)
    plt.gca().set_aspect("equal", adjustable="box")

    plt.title(f"{m}-QAM constellation")
    plt.xlabel("I")
    plt.ylabel("Q")

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_path = output_dir / f"qam_{m}.png"
        plt.savefig(output_path, dpi=300, bbox_inches="tight")

        print(f"Saved: {output_path}")

    if show:
        plt.show()
    else:
        plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot QAM constellations from txt file"
    )

    parser.add_argument(
        "filename",
        help="Path to txt file with constellation points"
    )

    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory where PNG images will be saved"
    )

    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not show plots, only save them"
    )

    args = parser.parse_args()

    constellations = read_constellations(args.filename)

    for m, points in constellations:
        plot_constellation(
            m,
            points,
            output_dir=args.output_dir,
            show=not args.no_show
        )


if __name__ == "__main__":
    main()