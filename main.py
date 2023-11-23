from typing import List
import sys
import os
from plot_scheme import plot_dynamic_of_process


if __name__ == "__main__":
    assert len(sys.argv) == 3, "нужно передать путь к файлу с результатом вычисления и название картинки"
    file_abspath = os.path.abspath(sys.argv[1])
    img_name = os.path.abspath(sys.argv[2])
    with open(file_abspath) as file:
        params = file.readline().split(' ')
        assert len(params) == 9
        I, K = list(map(int, params[:2]))
        T, l, k, c, alpha, beta, R = list(map(float, params[2:]))

        w_i_n: List[List[float]] = []
        for row_num in range(K + 1):
            row = file.readline().split(' ')
            w_i_n.append(list(map(float, row[:len(row) - 1])))

        assert len(w_i_n) == K + 1
        assert len(w_i_n[0]) == I + 1

        fig = plot_dynamic_of_process(w_i_n, I, K, l, T, True)
        fig.savefig(f"{img_name}.png")
