#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "mpi.h"

using namespace std;

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

int main(int argc, char* argv[]) {
    int ProcRank, ProcNum;

    double start, end;

    double al, be, ks, et;

    // константы
    int I           = 200;
    int K           = 200;
    int INTEGRAL    = 144;
    double T        = 150;
    double l        = 4;
    double k        = 0.01;
    double c        = 1.65;
    double alpha    = 0.005;
    double beta     = 0.25;
    double R        = 3;

    double p_i0;
    double q_i0;
    double ksi_i0;
    double eta_i0;

    const double h_z = l / I;
    const double h_t = T / K;

    // константы для основного уравнения
    const double gamma = k * h_t / (c * pow(h_z, 2));
    const double mu = 2 * alpha * h_t / (c * R);
    const double nu = 2 * h_t * INTEGRAL / (c * pow(R, 2));

    // константы для краевых условий
    const double var_theta = alpha * h_z / k;
    const double rho = c * pow(h_z, 2) / (2 * k * h_t);
    const double sigma = alpha * pow(h_z, 2) / (k * R);
    const double tau = pow(h_z, 2) * INTEGRAL / (k * pow(R, 2));

    int i_0 = I / 2;
    double **w_i_n_left, **w_i_n_right, **w_i_n;

    MPI_Status Status;

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

    double tbeg = MPI_Wtime();

    if (ProcRank == 0) {
        // cout << "First thread started work!" << endl;
        
        w_i_n_right = alloc_2d_double(K + 1, I + 1);
        fill_2d_array(w_i_n_right, K + 1, I + 1);

        // правые прогоночные коэффициенты
        double* p_i = new double[I];
        for (size_t i = 0; i < I; ++i) p_i[i] = 0.0;

        double** q_i_n = new double* [K + 1];
        for (size_t i = 0; i < K + 1; ++i) {
            q_i_n[i] = new double[I];
            for (size_t j = 0; j < I; ++j) q_i_n[i][j] = 0.0;
        }

        p_i[0] = 1 / (1 + var_theta + rho + sigma);
        for (int i = 1; i <= i_0; ++i) p_i[i] = gamma / (1 + mu + gamma * (2 - p_i[i - 1]));

        for (int n = 1; n < K + 1; ++n) {
            // cout << "proc #" << ProcRank << " n=" << n << endl;
            q_i_n[n][0] = (rho * w_i_n_right[n - 1][0] + tau) / (1 + var_theta + rho + sigma);
            // прямой правый ход
            for (int i = 1; i <= i_0; ++i) {
                double numerator = w_i_n_right[n - 1][i] + gamma * q_i_n[n][i - 1] + nu * exp(-beta * i * h_z);
                double denominator = 1 + mu + gamma * (2 - p_i[i - 1]);
                q_i_n[n][i] = numerator / denominator;
            }
            
            p_i0 = p_i[i_0];
            q_i0 = q_i_n[n][i_0];

            MPI_Send(&p_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&ksi_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
            MPI_Send(&q_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&eta_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);

            w_i_n_right[n][i_0] = (p_i[i_0] * eta_i0 + q_i_n[n][i_0]) / (1 - p_i[i_0] * ksi_i0);
            
            // обратный правый ход
            for (int i = i_0 - 1; i > -1; --i) w_i_n_right[n][i] = p_i[i] * w_i_n_right[n][i + 1] + q_i_n[n][i];
        }

        delete[] p_i;
        // cout << "p_i deleted!" << endl;
        for (size_t i = 0; i < K + 1; ++i) delete[] q_i_n[i];
        delete[] q_i_n;
        // cout << "q_i deleted!" << endl;

        w_i_n_left = alloc_2d_double(K + 1, I + 1);
        MPI_Recv(&(w_i_n_left[0][0]), (K + 1) * (I + 1), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
        // cout << "Left part of w_i_n recieved in first thread!" << endl;

        w_i_n = alloc_2d_double(K + 1, I + 1);

        for (int n = 1; n < K + 1; ++n) {
            for (int i = 0; i < i_0; ++i) w_i_n[n][i] = w_i_n_right[n][i];
            for (int i = i_0; i < I + 1; ++i) w_i_n[n][i] = w_i_n_left[n][i];
        }

        // cout << "Local w_i_n of first thread deleted!" << endl;
        double elapsedTime = MPI_Wtime() - tbeg;
        double totalTime;
        MPI_Reduce(&elapsedTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        cout << "MPI time: " << totalTime << endl;

        ofstream fout;
        fout.open("mpi_result.txt");
        if (!fout.is_open()) cout << "Ошибка при создании/открытии файла" << endl;
        else {
            fout << I << ' ' << K << ' ' << T << ' ' << l << ' ' << k << ' ' << c << ' ';
            fout << alpha << ' ' << beta << ' ' << R << endl;

            for (int i = 0; i < K + 1; ++i) {
                for (int j = 0; j < I + 1; ++j)
                    fout << w_i_n[i][j] << ' ';
                fout << endl;
            }
        }
        fout.close();

        free(w_i_n_right[0]);
        free(w_i_n_right);

        free(w_i_n_left[0]);
        free(w_i_n_left);

        free(w_i_n[0]);
        free(w_i_n);
    }
    
    if (ProcRank == 1) {
        // cout << "Second thread started work!" << endl;

        w_i_n_left = alloc_2d_double(K + 1, I + 1);
        fill_2d_array(w_i_n_left, K + 1, I + 1);

        // левые прогоночные коэффициенты
        double* ksi_i = new double[I + 1];
        for (size_t i = 0; i < I + 1; ++i) ksi_i[i] = 0.0;

        double** eta_i_n = new double* [K + 1];
        for (size_t i = 0; i < K + 1; ++i) {
            eta_i_n[i] = new double[I + 1];
            for (size_t j = 0; j < I + 1; ++j) eta_i_n[i][j] = 0.0;
        }

        ksi_i[I] = 1 / (1 + var_theta + rho + sigma);
        for (int i = I - 1; i >= i_0; --i) ksi_i[i] = gamma / (1 + mu + gamma * (2 - ksi_i[i + 1]));

        for (int n = 1; n < K + 1; ++n) {
            // cout << "proc #" << ProcRank << " n=" << n << endl;
            eta_i_n[n][I] = (rho * w_i_n_left[n - 1][I] + tau * exp(-beta * I * h_z)) / (1 + var_theta + rho + sigma);
            
            // прямой левый ход
            for (int i = I - 1; i >= i_0; --i) {
                double numerator = w_i_n_left[n - 1][i] + gamma * eta_i_n[n][i + 1] + nu * exp(-beta * i * h_z);
                double denominator = 1 + mu + gamma * (2 - ksi_i[i + 1]);
                eta_i_n[n][i] = numerator / denominator;
            }

            ksi_i0 = ksi_i[i_0 + 1];
            eta_i0 = eta_i_n[n][i_0 + 1];

            MPI_Send(&ksi_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&p_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
            MPI_Send(&eta_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&q_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);

            w_i_n_left[n][i_0] = (p_i0 * eta_i0 + q_i0) / (1 - p_i0 * ksi_i0);
            
            // обратный левый ход
            for (int i = i_0; i < I; ++i) w_i_n_left[n][i + 1] = ksi_i[i + 1] * w_i_n_left[n][i] + eta_i_n[n][i + 1];
        }
        delete[] ksi_i;
        // cout << "ksi_i deleted!" << endl;
        for (size_t i = 0; i < K + 1; ++i) delete[] eta_i_n[i];
        delete[] eta_i_n;
        // cout << "eta_i_n deleted!" << endl;

        MPI_Send(&(w_i_n_left[0][0]), (K + 1) * (I + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        // cout << "Left part of w_i_n sended from second thread!" << endl;

        free(w_i_n_left[0]);
        free(w_i_n_left);
        double elapsedTime = MPI_Wtime() - tbeg;
        double totalTime;
        MPI_Reduce(&elapsedTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}