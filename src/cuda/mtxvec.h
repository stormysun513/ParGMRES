#ifndef __MTXVEC_H__
#define __MTXVEC_H__

#include <vector>
#include <map>
#include <cmath>

class Vector
{
private:
    float* data;
    size_t len;
public:
    Vector();
    Vector(size_t size);

    void set(size_t idx, double val);
    double get(size_t idx) const;
    size_t size() const;
    float* getData() const;
};

inline void Vector::set(size_t idx, double val) { data[idx] = val; }
inline double Vector::get(size_t idx) const { return data[idx]; }
inline size_t Vector::size() const { return len; }
inline float* Vector::getData() const { return data; }

class Matrix
{
private:
    float* data;

    size_t n_rows;
    size_t n_cols;
public:
    Matrix();
    Matrix(size_t m, size_t n);
    Matrix(size_t m, size_t n, float* data_);

    size_t nRows() const;
    size_t nCols() const;
    void set(size_t row_idx, size_t col_idx, double val);
    double get(size_t row_idx, size_t col_idx) const;
};


inline size_t Matrix::nRows() const { return n_rows; }
inline size_t Matrix::nCols() const { return n_cols; }

inline void Matrix::set(size_t row_idx, size_t col_idx, double val) {
    data[row_idx*n_cols + col_idx] = val;
}

inline double Matrix::get(size_t row_idx, size_t col_idx) const {
    return data[row_idx*n_cols + col_idx];
}

#endif  // MTXVEC_H
