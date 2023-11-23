#include "utils.h"

double** alloc_2d_double(int rows, int cols) {
    double* data = (double *)malloc(rows * cols * sizeof(double));
    double** array = (double **)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; ++i)
        array[i] = &(data[cols * i]);

    return array;
}

void fill_2d_array(double** arr, int rows, int cols) {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) arr[i][j] = 0.0;
    }
}

void print_first(double** w_i_n, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << w_i_n[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void print_last(double** w_i_n, int rows, int cols, int n) {
    for (int i = rows - 1; i > rows - n; --i) {
        for (int j = cols - 1; j > cols - n; --j) {
            cout << w_i_n[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}

void write_to_file(ofstream &fout, double** w_i_n, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j)
            fout << w_i_n[i][j] << ' ';
        fout << endl;
    }
}