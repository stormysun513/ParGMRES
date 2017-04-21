#ifndef __MTXVEC_H__
#define __MTXVEC_H__

#include <vector>
#include <tuple>
#include <map>
#include <cmath>

class Matrix;

class Vector
{
private:
    std::vector<double> data;
public:
    Vector();
    Vector(size_t size);

    void resize(size_t size);
    void set(size_t idx, double val);
    double get(size_t idx) const;
    size_t size() const;
    void copy(const Vector& other);

    double norm2() const;
    double dotV(const Vector& other) const;
    Matrix crossV(const Vector& other) const;
    Vector& normalize();
    Vector& inverse();

    Vector add(const Vector& other) const;
    Vector& iadd(const Vector& other); // same as add, but inplace
    Vector sub(const Vector& other) const;
    Vector& isub(const Vector& other); // same as sub, but inplace
    Vector mulS(double scaler) const;
    Vector& imulS(double scaler); // same as mulS, but inplace
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
    Vector getRow(size_t row_idx) const;
    Vector getCol(size_t col_idx) const;
    Vector mul(const Vector& vec) const;
    Vector mulPartial(const Vector& vec, size_t n_cols_) const;
    Vector mulPartialT(const Vector& vec, size_t n_rows_) const;

    Matrix& isub(const Matrix& other);
    Matrix& iadd(const Matrix& other);
    Matrix& imulS(double scaler); // same as mulS, but inplace
    Matrix& iRowMulS(const Vector& other);
    Matrix& iRowDivS(const Vector& other);

    Matrix covariance() const;
    Matrix transpose() const;
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

// CSC format
class SparseMatrix
{
private:
    std::vector<double> data;

    // row indices
    std::vector<size_t> indices;

    // col start
    std::vector<size_t> indptr;

    size_t n_rows;
    size_t n_cols;

    bool posInData(size_t i, size_t j, size_t& ret) const;
    void construct(std::vector<std::tuple<double, size_t, size_t>> raw_data,
                   size_t n_rows, size_t n_cols);
public:
    SparseMatrix(std::vector<std::tuple<double, size_t, size_t>> raw_data,
                 size_t n_rows, size_t n_cols);

    SparseMatrix(Matrix dense);

    size_t nnz() const;
    size_t nRows() const;
    size_t nCols() const;
    double get(size_t row_idx, size_t col_idx) const;
    bool set(size_t row_idx, size_t col_idx, double val);

    Vector getCol(size_t col_idx) const;
    Vector mul(const Vector& vec) const;
    Vector mulPartial(const Vector& vec, size_t n_cols_) const;
};

inline double SparseMatrix::get(size_t i, size_t j) const {
    size_t idx;
    if (posInData(i, j, idx))
        return data[idx];
    else
        return 0.0;
}

inline bool SparseMatrix::set(size_t i, size_t j, double val) {
    size_t idx;
    if (posInData(i, j, idx)) {
        data[idx] = val;
        return true;
    } else {
        return false;
    }
}

inline size_t SparseMatrix::nnz() const { return data.size(); }
inline size_t SparseMatrix::nRows() const { return n_rows; }
inline size_t SparseMatrix::nCols() const { return n_cols; }

#endif  // MTXVEC_H
