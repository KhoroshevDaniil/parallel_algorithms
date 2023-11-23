#pragma once
#include <iostream>
#include <fstream>

using namespace std;

double** alloc_2d_double(int rows, int cols);
void fill_2d_array(double** arr, int rows, int cols);
void print_first(double** w_i_n, int n);
void print_last(double** w_i_n, int rows, int cols, int n);
void write_to_file(ofstream &fout, double** w_i_n, int rows, int cols);