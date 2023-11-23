from typing import List
import numpy as np
from matplotlib import pyplot as plt


def to_fixed(numObj, digits=0):
    return f"{numObj:.{digits}f}"


def plot_dynamic_of_process(w: List[List[float]], I: int, K: int, l: float, T: float, show=False):
    w = np.array(w)
    x1 = np.linspace(0, l, I + 1)
    x2 = np.linspace(0, T, K + 1)
    z_arr = np.linspace(0, I, 8, dtype=int)
    t_arr = np.linspace(0, K, 8, dtype=int)

    h_z = l / I
    h_t = T / K

    fig, (ax1, ax2) = plt.subplots(1, 2)

    for t in t_arr:
        ax1.plot(x1, w[t], label=f"t = {to_fixed(t * h_t, 1)}")

    for z in z_arr:
        ax2.plot(x2, w[:, z], label=f"z = {to_fixed(z * h_z, 1)}")

    fig.suptitle(f"Неявная модифицированная схема\n{I=}, {K=}")

    ax1.set_xlabel('Точка на стержне, см')
    ax1.set_ylabel('Температура, град')
    ax1.set_title('Распределение температуры по всему стержню')

    ax2.set_xlabel('Время, с')
    ax2.set_ylabel('Температура, град')
    ax2.set_title('Динамика температуры в точках стержня')

    for ax in ax1, ax2:
        ax.grid(True)
        ax.legend()

    fig.set_size_inches(10, 6)
    plt.tight_layout()

    if show:
        plt.show()

    return fig
