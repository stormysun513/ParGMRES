#ifndef __MTXVEC_H__
#define __MTXVEC_H__

#include <vector>
#include <cmath>

class Vector
{
private:
    std::vector<double> data;
public:
    Vector(size_t size);

    void set(size_t idx, double val);
    double get(size_t idx) const;
    size_t size() const;

    double norm2() const;
    double dotV(const Vector& other) const;

    Vector add(const Vector& other) const;
    Vector iadd(const Vector& other); // same as add, but inplace
    Vector sub(const Vector& other) const;
    Vector isub(const Vector& other); // same as sub, but inplace
    Vector mulS(double scaler) const;
    Vector imulS(double scaler); // same as mulS, but inplace
};

class Matrix
{
private:
    std::vector<std::vector<double>> data;
    size_t n_rows;
    size_t n_cols;
public:
    Matrix();
    Matrix(size_t m, size_t n);

    void resize(size_t m, size_t n);
    size_t nRows() const;
    size_t nCols() const;
    void set(size_t row_idx, size_t col_idx, double val);
    double get(size_t row_idx, size_t col_idx) const;

    void setRow(size_t row_idx, const Vector& vec);
    void setCol(size_t col_idx, const Vector& vec);
    Vector getCol(size_t col_idx) const;
    Vector mul(const Vector& vec) const;
    Vector mulPartial(const Vector& vec, size_t n_cols_) const;
};


inline void Vector::set(size_t idx, double val) { data[idx] = val; }
inline double Vector::get(size_t idx) const { return data[idx]; }
inline size_t Vector::size() const { return data.size(); }
inline size_t Matrix::nRows() const { return n_rows; }
inline size_t Matrix::nCols() const { return n_cols; }

inline void Matrix::set(size_t row_idx, size_t col_idx, double val) {
    data[row_idx][col_idx] = val;
}

inline double Matrix::get(size_t row_idx, size_t col_idx) const {
    return data[row_idx][col_idx];
}

#endif  // MTXVEC_H
