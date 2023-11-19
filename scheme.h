#pragma once
#include <math.h>
#include <vector>
#include <time.h>
#include <iostream>

using namespace std;

const int INTEGRAL = 144;

// w_i_n - returning 2D array, K + 1 x I + 1
//vector<vector<double>> calculate_implicit_modified_scheme(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R);
double** right_run_through(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R);
double** left_run_through(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R);
double** counter_run_through(const int& K, const int& I, const double& l, const double& T, const double& k, const double& c, const double& alpha, const double& beta, const double& R);