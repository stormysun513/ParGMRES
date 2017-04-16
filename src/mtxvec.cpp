#include "mtxvec.h"

#include <vector>
#include <cmath>
#include <cassert>

/*
 * Vector member funciton
 */
Vector::Vector(size_t size) {
    data.resize(size);
}

double Vector::norm2() const {
    double res = .0f;

    for (double num: data) {
        res += num * num;
    }
    return sqrt(res);
}

double Vector::dotV(const Vector& other) const {
    double sum = .0f;

    assert(data.size() == other.size());

    for (size_t i = 0; i < data.size(); i++) {
        sum += data[i] * other.get(i);
    }

    return sum;
}

Vector Vector::add(const Vector& other) const {
    assert(data.size() == other.size());

    Vector res(data.size());

    for (size_t i = 0; i < data.size(); i++) {
        res.set(i, data[i] + other.get(i));
    }

    return res;
}

// same as add, but inplace
Vector Vector::iadd(const Vector& other) {
    assert(data.size() == other.size());

    for (size_t i = 0; i < data.size(); i++) {
        data[i] += other.get(i);
    }

    return *this;
}

Vector Vector::sub(const Vector& other) const {
    assert(data.size() == other.size());

    Vector res(data.size());

    for (size_t i = 0; i < data.size(); i++) {
        res.set(i, data[i] - other.get(i));
    }

    return res;
}

// same as sub, but inplace
Vector Vector::isub(const Vector& other) {
    assert(data.size() == other.size());

    for (size_t i = 0; i < data.size(); i++) {
        data[i] -= other.get(i);
    }

    return *this;
}

Vector Vector::mulS(double scaler) const {
    Vector res(data.size());

    for (size_t i = 0; i < data.size(); i++) {
        res.set(i, data[i] * scaler);
    }

    return res;
}

// same as mulS, but inplace
Vector Vector::imulS(double scaler) {
    for (size_t i = 0; i < data.size(); i++) {
        data[i] *= scaler;
    }

    return *this;
}

/*
 * Matrix member funciton
 */
Matrix::Matrix():n_rows(0),n_cols(0){}

Matrix::Matrix(size_t m, size_t n) {
    data.resize(m);
    n_rows = m;
    n_cols = n;
    for (int i = 0; i < m; i++){
        data[i] = std::vector<double>(n);
    }
}

void Matrix::resize(size_t m, size_t n){
    data.resize(m);
    n_rows = m;
    n_cols = n;
    for(int i = 0; i < m; i++) {
        data[i] = std::vector<double>(n);
    }
}

void Matrix::setRow(size_t row_idx, const Vector& vec) {
    size_t len = vec.size();

    assert(0 <= row_idx);
    assert(row_idx < n_rows);
    assert(len == n_cols);

    for (size_t i = 0; i < len; ++i) {
        data[row_idx][i] = vec.get(i);
    }
}

void Matrix::setCol(size_t col_idx, const Vector& vec) {
    size_t len = vec.size();

    assert(0 <= col_idx);
    assert(col_idx < n_cols);
    assert(len == n_rows);

    for (size_t i = 0; i < len; ++i) {
        data[i][col_idx] = vec.get(i);
    }
}

Vector Matrix::getCol(size_t col_idx) const {
    Vector col(n_rows);

    for (int i = 0; i < n_rows; ++i) {
        col.set(i, data[i][col_idx]);
    }

    return col;
}

Vector Matrix::mul(const Vector& vec) const {
    Vector ret(n_rows);

    for (size_t i = 0; i < n_rows; ++i) {
        ret.set(i, 0.0);
        for (size_t j = 0; j < n_cols; ++j) {
            ret.set(i, ret.get(i) + data[i][j] * vec.get(j));
        }
    }

    return ret;
}

Vector Matrix::mulPartial(const Vector& vec, size_t n_cols_) const {
    assert(n_cols >= vec.size());
    Vector ret(n_rows);

    for (size_t i = 0; i < n_rows; ++i) {
        ret.set(i, 0.0);
        for (size_t j = 0; j < n_cols_; ++j) {
            ret.set(i, ret.get(i) + data[i][j] * vec.get(j));
        }
    }

    return ret;
}
