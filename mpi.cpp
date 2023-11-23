#include "mpi.h"
#include <time.h>
#include <math.h>
#include "utils.h"


int main(int argc, char* argv[]) {
    int ProcRank, ProcNum;

    // регулируемые параметры
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

    // прогоночные коэффициенты, которыми процессы обмениваются друг с другом
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

    // точка встречной прогонки (по умолчанию - серидинка)
    int i_0 = I / 2;

    // матрицы левой и правой прогонок 
    double **w_i_n_left, **w_i_n_right;

    MPI_Status Status;

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

    double tbeg = MPI_Wtime();

    if (ProcRank == 0) {
        // создаём матрицу правой прогонки и инициалищируем её нулями
        w_i_n_right = alloc_2d_double(K + 1, I + 1);
        fill_2d_array(w_i_n_right, K + 1, I + 1);

        // то же самое делаем для прогоночных коэффициентов
        double* p_i = new double[I];
        for (size_t i = 0; i < I; ++i) p_i[i] = 0.0;

        double** q_i_n = alloc_2d_double(K + 1, I);
        fill_2d_array(q_i_n, K + 1, I);

        // начальное значение прогоночного коэффициента взято из левого краевого условия
        p_i[0] = 1 / (1 + var_theta + rho + sigma);
        
        // далее вычисляем прогоночные коэффициенты, которые не зависят от слоя по времени "n" (одинаковы для всех "n")
        for (int i = 1; i <= i_0; ++i) p_i[i] = gamma / (1 + mu + gamma * (2 - p_i[i - 1]));

        for (int n = 1; n < K + 1; ++n) {
            // начальное значение взято из левого краевого условия
            q_i_n[n][0] = (rho * w_i_n_right[n - 1][0] + tau) / (1 + var_theta + rho + sigma);
           
            // прямой правый ход
            for (int i = 1; i <= i_0; ++i) {
                double numerator = w_i_n_right[n - 1][i] + gamma * q_i_n[n][i - 1] + nu * exp(-beta * i * h_z);
                double denominator = 1 + mu + gamma * (2 - p_i[i - 1]);
                q_i_n[n][i] = numerator / denominator;
            }
            
            p_i0 = p_i[i_0];
            q_i0 = q_i_n[n][i_0];

            // обмен прогоночными коэффициентами
            MPI_Send(&p_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&ksi_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);
            MPI_Send(&q_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&eta_i0, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);

            w_i_n_right[n][i_0] = (p_i[i_0] * eta_i0 + q_i_n[n][i_0]) / (1 - p_i[i_0] * ksi_i0);
            
            // обратный правый ход
            for (int i = i_0 - 1; i > -1; --i) w_i_n_right[n][i] = p_i[i] * w_i_n_right[n][i + 1] + q_i_n[n][i];
        }

        delete[] p_i;
        free(q_i_n[0]);
        free(q_i_n);

        // от второй задачи получаем матрицу левой прогонки
        w_i_n_left = alloc_2d_double(K + 1, I + 1);
        MPI_Recv(&(w_i_n_left[0][0]), (K + 1) * (I + 1), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);

        // дополняем матрицу правой прогонки значениями матрицы левой прогонки
        for (int n = 1; n < K + 1; ++n) {
            for (int i = i_0; i < I + 1; ++i) w_i_n_right[n][i] = w_i_n_left[n][i];
        }

        double elapsedTime = MPI_Wtime() - tbeg;
        double totalTime;
        MPI_Reduce(&elapsedTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        cout << "MPI time: " << totalTime << endl;

        // записываем матрицу в файл        
        ofstream fout;
        fout.open("mpi_result.txt");
        if (!fout.is_open()) cout << "Ошибка при создании/открытии файла" << endl;
        else {
            fout << I << ' ' << K << ' ' << T << ' ' << l << ' ' << k << ' ' << c << ' ';
            fout << alpha << ' ' << beta << ' ' << R << endl;

            write_to_file(fout, w_i_n_right, K + 1, I + 1);
        }
        fout.close();

        free(w_i_n_right[0]);
        free(w_i_n_right);

        free(w_i_n_left[0]);
        free(w_i_n_left);
    }
    
    if (ProcRank == 1) {
        // создаём матрицу левой прогонки и инициалищируем её нулями
        w_i_n_left = alloc_2d_double(K + 1, I + 1);
        fill_2d_array(w_i_n_left, K + 1, I + 1);

        // то же самое делаем для прогоночных коэффициентов
        double* ksi_i = new double[I + 1];
        for (size_t i = 0; i < I + 1; ++i) ksi_i[i] = 0.0;

        double** eta_i_n = alloc_2d_double(K + 1, I + 1);
        fill_2d_array(eta_i_n, K + 1, I + 1);

        // начальное значение прогоночного коэффициента взято из правого краевого условия
        ksi_i[I] = 1 / (1 + var_theta + rho + sigma);
        // далее вычисляем прогоночные коэффициенты, которые не зависят от слоя по времени "n" (одинаковы для всех "n") 
        for (int i = I - 1; i >= i_0; --i) ksi_i[i] = gamma / (1 + mu + gamma * (2 - ksi_i[i + 1]));

        for (int n = 1; n < K + 1; ++n) {
            eta_i_n[n][I] = (rho * w_i_n_left[n - 1][I] + tau * exp(-beta * I * h_z)) / (1 + var_theta + rho + sigma);
            
            // прямой левый ход
            for (int i = I - 1; i >= i_0; --i) {
                double numerator = w_i_n_left[n - 1][i] + gamma * eta_i_n[n][i + 1] + nu * exp(-beta * i * h_z);
                double denominator = 1 + mu + gamma * (2 - ksi_i[i + 1]);
                eta_i_n[n][i] = numerator / denominator;
            }

            ksi_i0 = ksi_i[i_0 + 1];
            eta_i0 = eta_i_n[n][i_0 + 1];

            // обмен прогоночными коэффициентами
            MPI_Send(&ksi_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&p_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
            MPI_Send(&eta_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Recv(&q_i0, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);

            w_i_n_left[n][i_0] = (p_i0 * eta_i0 + q_i0) / (1 - p_i0 * ksi_i0);
            
            // обратный левый ход
            for (int i = i_0; i < I; ++i) w_i_n_left[n][i + 1] = ksi_i[i + 1] * w_i_n_left[n][i] + eta_i_n[n][i + 1];
        }
        delete[] ksi_i;
        free(eta_i_n[0]);
        free(eta_i_n);

        // отправляем первой задаче вычисленные значения матрицы левой прогонки
        MPI_Send(&(w_i_n_left[0][0]), (K + 1) * (I + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

        double elapsedTime = MPI_Wtime() - tbeg;
        double totalTime;
        MPI_Reduce(&elapsedTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        free(w_i_n_left[0]);
        free(w_i_n_left);
    }

    MPI_Finalize();
    return 0;
}