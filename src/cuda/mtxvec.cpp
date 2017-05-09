#include "mtxvec.h"

#include <vector>

#include <iostream>
using namespace std;

/*
 * Vector member funciton
 */
Vector::Vector(){}

Vector::Vector(size_t size) {
    len = size;
    data = (float *) calloc(size, sizeof(float));
}

/*
 * Matrix member funciton
 */
Matrix::Matrix():n_rows(0),n_cols(0){}

Matrix::Matrix(size_t m, size_t n) {
    n_rows = m;
    n_cols = n;
    data = (float *)calloc(m * n, sizeof(float));
}

Matrix::Matrix(size_t m, size_t n, float* data_) {
    n_rows = m;
    n_cols = n;
    data = data_;
}
