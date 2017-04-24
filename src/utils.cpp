#include <iostream>
#include <random>
#include <tuple>
#include <cassert>

#include "utils.h"

void
printVector(const Vector& vec){
    size_t length = vec.size();

    char buf[1024];
    buf[0] = '\0';

    for (size_t j = 0; j < length; j++) {
        if (j == 0)
            sprintf(buf, "%s[", buf);

        sprintf(buf, "%s%.4f", buf, vec.get(j));
        if (j == length-1)
            sprintf(buf, "%s]", buf);
        else
            sprintf(buf, "%s, ", buf);
    }
    sprintf(buf, "%s\n", buf);
    std::cout << buf;
}

void
printMatrix(const Matrix& mat,
            size_t row_start, size_t row_end,
            size_t col_start, size_t col_end) {
    size_t m = mat.nRows();
    size_t n = mat.nCols();

    assert(row_start >= 0);
    assert(col_start >= 0);
    assert(row_end <= m);
    assert(col_end <= n);

    std::cout << "[";
    for (size_t i = row_start; i < row_end; i++) {
        for (size_t j = col_start; j < col_end; j++) {
            if (j == col_start)
                std::cout << "[";

            std::cout << mat.get(i, j);
            if (j == col_end-1)
                std::cout << "],";
            else
                std::cout << ", ";
        }
        std::cout << "\n";
    }
    std::cout << "]\n";
}

Vector randUniformVector(size_t n){
    Vector vec(n);
    std::random_device rd;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    for(size_t i = 0; i < n; i++){
        vec.set(i, distribution(rd));
    }
    return vec;
}

Vector randUnitUniformVector(size_t n){
    Vector vec(n);
    std::random_device rd;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    for(size_t i = 0; i < n; i++){
        vec.set(i, distribution(rd));
    }
    return vec.normalize();
}

CSRMatrix
identitySparseMatrix(size_t size) {

    std::vector<std::tuple<double, size_t, size_t>> raw_data;
    for (size_t i = 0; i < size; ++i) {
        raw_data.push_back(std::make_tuple(1, i, i));
    }

    CSRMatrix mat_csr(raw_data, size, size);
    return mat_csr;
}

Matrix
identityMatrix(size_t size) {
    Matrix mat(size, size);

    for (size_t i = 0; i < size; i++) {
        mat.set(i, i, 1.0);
    }
    return mat;
}
