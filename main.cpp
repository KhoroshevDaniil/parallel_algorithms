#include "scheme.h"
#include "utils.h"


int main(int argc, char* argv[]) {
    int I           = 200;
    int K           = 200;
    double T        = 150;
    double l        = 4;
    double k        = 0.01;
    double c        = 1.65;
    double alpha    = 0.005;
    double beta     = 0.25;
    double R        = 3;

    //vector<vector<double>> w_i_n = calculate_implicit_modified_scheme(K, I, l, T, k, c, alpha, beta, R);
    double** w_i_n_right    = right_run_through(K, I, l, T, k, c, alpha, beta, R);
    double** w_i_n_left     = left_run_through(K, I, l, T, k, c, alpha, beta, R);
    double** w_i_n_counter  = counter_run_through(K, I, l, T, k, c, alpha, beta, R);

    int rows_count = K + 1;
    int cols_count = I + 1;


    cout << "First 5 from RIGHT run-through algorithm" << endl;
    print_first(w_i_n_right, 5);

    cout << "First 5 from LEFT run-through algorithm" << endl;
    print_first(w_i_n_left, 5);

    cout << "First 5 from COUNTER run-through algorithm" << endl;
    print_first(w_i_n_counter, 5);


    cout << "Last 5 from RIGHT run-through algorithm" << endl;
    print_last(w_i_n_right, rows_count, cols_count, 5);

    cout << "Last 5 from LEFT run-through algorithm" << endl;
    print_last(w_i_n_left, rows_count, cols_count, 5);

    cout << "Last 5 from COUNTER run-through algorithm" << endl;
    print_last(w_i_n_counter, rows_count, cols_count, 5);


    double error_right_left = 0.0;
    double error_right_counter = 0.0;
    for(int i = 0; i < rows_count; ++i) {
        for (int j = 0; j < cols_count; ++j) {
            error_right_left    += w_i_n_right[i][j] - w_i_n_left[i][j];
            error_right_counter += w_i_n_right[i][j] - w_i_n_counter[i][j];
        }
    }

    cout << "Error between RIGHT and LEFT:\t\t" << error_right_left << endl;
    cout << "Error between RIGHT and COUNTER:\t" << error_right_counter << endl;

    ofstream fout;

    fout.open("cpu_result.txt");
    if (!fout.is_open()) cout << "Ошибка при создании/открытии файла" << endl;
    else {
        fout << I << ' ' << K << ' ' << T << ' ' << l << ' ' << k << ' ' << c << ' ';
        fout << alpha << ' ' << beta << ' ' << R << endl;

        write_to_file(fout, w_i_n_counter, rows_count, cols_count);
    }
    fout.close();

    for (size_t i = 0; i < K + 1; ++i) delete[] w_i_n_right[i];
    delete[] w_i_n_right;

    for (size_t i = 0; i < K + 1; ++i) delete[] w_i_n_left[i];
    delete[] w_i_n_left;

    for (size_t i = 0; i < K + 1; ++i) delete[] w_i_n_counter[i];
    delete[] w_i_n_counter;

    return 0;
}