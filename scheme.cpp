#include "scheme.h"

double** right_run_through(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R) {
    double start, end;

    //std::vector<double> p_i(I);
    double* p_i = new double[I];
    for (size_t i = 0; i < I; ++i) p_i[i] = 0.0;

    //std::vector<std::vector<double>> q_i_n(K + 1, std::vector<double>(I));
    double** q_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        q_i_n[i] = new double[I];
        for (size_t j = 0; j < I; ++j) q_i_n[i][j] = 0.0;
    }

    //std::vector<std::vector<double>> w_i_n(K + 1, std::vector<double>(I + 1));
    double** w_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        w_i_n[i] = new double[I + 1];
        for (size_t j = 0; j < I + 1; ++j) w_i_n[i][j] = 0.0;
    }

    start = clock();

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

    // p_0 got from left bound condition
    p_i[0] = 1 / (1 + var_theta + rho + sigma);

    // calculating p_i coefficients, that the same for each n layer
    for (int i = 1; i < I; ++i) p_i[i] = gamma / (1 + mu + gamma * (2 - p_i[i - 1]));

    for (int n = 1; n < K + 1; ++n) {
        // q_0_n got from left bound 
        q_i_n[n][0] = (rho * w_i_n[n - 1][0] + tau) / (1 + var_theta + rho + sigma);

        for (int i = 1; i < I; ++i) {
            double numerator = w_i_n[n - 1][i] + gamma * q_i_n[n][i - 1] + nu * exp(-beta * i * h_z);
            double denominator = 1 + mu + gamma * (2 - p_i[i - 1]);
            q_i_n[n][i] = numerator / denominator;
        }

        // start calculate w_i_n from point i = I
        double w_I_n_numerator = rho * w_i_n[n - 1][I] + q_i_n[n][I - 1] + tau * exp(-beta * I * h_z);
        double w_I_n_denominator = 1 + var_theta + rho + sigma - p_i[I - 1];
        w_i_n[n][I] = w_I_n_numerator / w_I_n_denominator;

        for (int i = I - 1; i > -1; --i) w_i_n[n][i] = p_i[i] * w_i_n[n][i + 1] + q_i_n[n][i];
    }
    end = clock();

    std::cout << "Right execution time: " << (end - start) / CLOCKS_PER_SEC << endl << endl;

    delete[] p_i;
    for (size_t i = 0; i < K + 1; ++i) delete[] q_i_n[i];
    delete[] q_i_n;

    return w_i_n;
}

double** left_run_through(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R) {
    double start, end;

    double* ksi_i = new double[I + 1];
    for (size_t i = 0; i < I + 1; ++i) ksi_i[i] = 0.0;

    double** eta_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        eta_i_n[i] = new double[I + 1];
        for (size_t j = 0; j < I + 1; ++j) eta_i_n[i][j] = 0.0;
    }

    double** w_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        w_i_n[i] = new double[I + 1];
        for (size_t j = 0; j < I + 1; ++j) w_i_n[i][j] = 0.0;
    }

    start = clock();
    
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

    // ksi_I got from left bound condition
    ksi_i[I] = 1 / (1 + var_theta + rho + sigma);

    // calculating p_i coefficients, that the same for each n layer
    for (int i = I - 1; i > -1; --i) ksi_i[i] = gamma / (1 + mu + gamma * (2 - ksi_i[i + 1]));

    for (int n = 1; n < K + 1; ++n) {
        // q_I_n got from left bound 
        eta_i_n[n][I] = (rho * w_i_n[n - 1][I] + tau * exp(-beta * I * h_z)) / (1 + var_theta + rho + sigma);

        for (int i = I - 1; i > -1; --i) {
            double numerator = w_i_n[n - 1][i] + gamma * eta_i_n[n][i + 1] + nu * exp(-beta * i * h_z);
            double denominator = 1 + mu + gamma * (2 - ksi_i[i + 1]);
            eta_i_n[n][i] = numerator / denominator;
        }

        // start calculate w_i_n from point i = 0
        double w_0_n_numerator = rho * w_i_n[n - 1][0] + eta_i_n[n][1] + tau;
        double w_0_n_denominator = 1 + var_theta + rho + sigma - ksi_i[1];
        w_i_n[n][0] = w_0_n_numerator / w_0_n_denominator;

        for (int i = 0; i < I; ++i) w_i_n[n][i + 1] = ksi_i[i + 1] * w_i_n[n][i] + eta_i_n[n][i + 1];
    }
    end = clock();
    cout << "Left execution time: " << (end - start) / CLOCKS_PER_SEC << endl << endl;

    delete[] ksi_i;
    for (size_t i = 0; i < K + 1; ++i) delete[] eta_i_n[i];
    delete[] eta_i_n;

    return w_i_n;
}

double** counter_run_through(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R) {
    double start, end;
    int i_0 = I / 2;
    cout << "i_0 = " << i_0 << endl;

    double* p_i = new double[I];
    for (size_t i = 0; i < I; ++i) p_i[i] = 0.0;

    double** q_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        q_i_n[i] = new double[I];
        for (size_t j = 0; j < I; ++j) q_i_n[i][j] = 0.0;
    }

    double* ksi_i = new double[I + 1];
    for (size_t i = 0; i < I + 1; ++i) ksi_i[i] = 0.0;

    double** eta_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        eta_i_n[i] = new double[I + 1];
        for (size_t j = 0; j < I + 1; ++j) eta_i_n[i][j] = 0.0;
    }

    double** w_i_n = new double* [K + 1];
    for (size_t i = 0; i < K + 1; ++i) {
        w_i_n[i] = new double[I + 1];
        for (size_t j = 0; j < I + 1; ++j) w_i_n[i][j] = 0.0;
    }

    start = clock();
    
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

    p_i[0] = 1 / (1 + var_theta + rho + sigma);
    ksi_i[I] = 1 / (1 + var_theta + rho + sigma);

    for (int i = 1; i <= i_0; ++i) p_i[i] = gamma / (1 + mu + gamma * (2 - p_i[i - 1]));
    for (int i = I - 1; i >= i_0; --i) ksi_i[i] = gamma / (1 + mu + gamma * (2 - ksi_i[i + 1]));

    for (int n = 1; n < K + 1; ++n) {
        q_i_n[n][0] = (rho * w_i_n[n - 1][0] + tau) / (1 + var_theta + rho + sigma);
        eta_i_n[n][I] = (rho * w_i_n[n - 1][I] + tau * exp(-beta * I * h_z)) / (1 + var_theta + rho + sigma);
        
        // прямой правый ход
        for (int i = 1; i <= i_0; ++i) {
            double numerator = w_i_n[n - 1][i] + gamma * q_i_n[n][i - 1] + nu * exp(-beta * i * h_z);
            double denominator = 1 + mu + gamma * (2 - p_i[i - 1]);
            q_i_n[n][i] = numerator / denominator;
        }

        // прямой левый ход
        for (int i = I - 1; i >= i_0; --i) {
            double numerator = w_i_n[n - 1][i] + gamma * eta_i_n[n][i + 1] + nu * exp(-beta * i * h_z);
            double denominator = 1 + mu + gamma * (2 - ksi_i[i + 1]);
            eta_i_n[n][i] = numerator / denominator;
        }

        w_i_n[n][i_0] = (p_i[i_0] * eta_i_n[n][i_0 + 1] + q_i_n[n][i_0]) / (1 - p_i[i_0] * ksi_i[i_0 + 1]);
        
        // обратный правый ход
        for (int i = i_0 - 1; i > -1; --i) w_i_n[n][i] = p_i[i] * w_i_n[n][i + 1] + q_i_n[n][i];
        // обратный левый ход
        for (int i = i_0; i < I; ++i) w_i_n[n][i + 1] = ksi_i[i + 1] * w_i_n[n][i] + eta_i_n[n][i + 1];
    }
    end = clock();

    cout << "Counter execution time: " << (end - start) / CLOCKS_PER_SEC << endl << endl;

    delete[] p_i;
    for (size_t i = 0; i < K + 1; ++i) delete[] q_i_n[i];
    delete[] q_i_n;
    
    delete[] ksi_i;
    for (size_t i = 0; i < K + 1; ++i) delete[] eta_i_n[i];
    delete[] eta_i_n;

    return w_i_n;
}