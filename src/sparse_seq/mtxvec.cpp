#include "mtxvec.h"

#include <vector>
#include <cmath>
#include <cassert>
#include <tuple>
#include <algorithm>

#include <omp.h>

#include <iostream>
using namespace std;

#define OMP_N_BOUND     1000000
#define OMP_NN_BOUND    500

static void sortRawDataByRow(std::vector<std::tuple<double, size_t, size_t>>& raw_data);

/*
 * Vector member funciton
 */
Vector::Vector(){}

Vector::Vector(size_t size) {
    data.resize(size);
}

void Vector::clear(){
    std::fill(data.begin(), data.end(), 0);
}

void Vector::resize(size_t size){
    data.resize(size);
}

void Vector::copy(const Vector& other){

    size_t size = other.size();

    assert(size == this->size());

    for(size_t i = 0; i < size; i++){
        this->set(i, other.get(i));
    }
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
    size_t size = data.size();

    assert(size == other.size());

    for (size_t i = 0; i < size; i++) {
        sum += data[i] * other.get(i);
    }
    return sum;
}

Matrix Vector::crossV(const Vector& other) const {

    size_t m = this->size();
    size_t n = other.size();
    Matrix mat(m, n);

    for(size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            mat.set(i, j, this->get(i)*other.get(j));
        }
    }
    return mat;
}

Vector& Vector::normalize(){

    double norm = this->norm2();
    size_t size = data.size();

    assert(norm != 0);

    for(size_t i = 0; i < size; i++){
        data[i] /= norm;
    }
    return *this;
}

Vector& Vector::inverse(){

    size_t size = data.size();

    for(size_t i = 0; i < size; i++){
        data[i] = 1.0/data[i];
    }
    return *this;
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
Vector& Vector::iadd(const Vector& other) {
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
Vector& Vector::isub(const Vector& other) {
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
Vector& Vector::imulS(double scaler) {
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
    for (size_t i = 0; i < m; i++){
        data[i] = std::vector<double>(n);
    }
}

void Matrix::clear(){
    for(auto vec: data)
    std::fill(vec.begin(), vec.end(), 0);
}

void Matrix::resize(size_t m, size_t n){

    data.resize(m);

    if(m > n_rows){
        for(size_t i = n_rows; i < m; i++){
            data[i] = std::vector<double>(n);
        }
    }
    for(size_t i = 0; i < n_rows; i++){
        data[i].resize(n);
    }
    n_rows = m;
    n_cols = n;
}

void Matrix::setRow(size_t row_idx, const Vector& vec) {
    size_t len = vec.size();

    assert(row_idx < n_rows);
    assert(len == n_cols);

    for (size_t i = 0; i < len; ++i) {
        data[row_idx][i] = vec.get(i);
    }
}

void Matrix::setCol(size_t col_idx, const Vector& vec) {
    size_t len = vec.size();

    assert(col_idx < n_cols);
    assert(len == n_rows);

    for (size_t i = 0; i < len; ++i) {
        data[i][col_idx] = vec.get(i);
    }
}

Vector Matrix::getRow(size_t row_idx) const {
    Vector row(n_cols);

    for (size_t i = 0; i < n_cols; ++i) {
        row.set(i, data[row_idx][i]);
    }
    return row;
}

Vector Matrix::getCol(size_t col_idx) const {
    Vector col(n_rows);

    for (size_t i = 0; i < n_rows; ++i) {
        col.set(i, data[i][col_idx]);
    }
    return col;
}

Vector Matrix::mul(const Vector& vec) const {
    Vector ret(n_rows);

    for (size_t i = 0; i < n_rows; ++i) {
        double temp = .0f;
        for (size_t j = 0; j < n_cols; ++j) {
            temp += data[i][j] * vec.get(j);
        }
        ret.set(i, temp);
    }
    return ret;
}

Vector Matrix::mulPartial(const Vector& vec, size_t n_cols_) const {
    assert(n_cols_ == vec.size());
    Vector ret(n_rows);

    for (size_t i = 0; i < n_rows; ++i) {
        double temp = .0f;
        for (size_t j = 0; j < n_cols_; ++j) {
            temp += data[i][j] * vec.get(j);
        }
        ret.set(i, temp);
    }
    return ret;
}

Vector Matrix::mulPartialT(const Vector& vec, size_t n_rows_) const {
    assert(n_rows_ == vec.size());
    Vector ret(n_cols);

    for (size_t j = 0; j < n_cols; ++j) {
        double temp = .0f;
        for (size_t i = 0; i < n_rows_; ++i) {
            temp += data[i][j] * vec.get(i);
        }
        ret.set(j, temp);
    }
    return ret;
}

Matrix& Matrix::isub(const Matrix& other){

    size_t m = other.nRows();
    size_t n = other.nCols();

    assert(m == n_rows);
    assert(n == n_cols);

    for(size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            data[i][j] -= other.get(i, j);
        }
    }
    return *this;
}

Matrix& Matrix::iadd(const Matrix& other){

    size_t m = other.nRows();
    size_t n = other.nCols();

    assert(m == n_rows);
    assert(n == n_cols);

    for(size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            data[i][j] += other.get(i, j);
        }
    }
    return *this;
}

Matrix& Matrix::imulS(double scalar){

    for(size_t i = 0; i < n_rows; i++){
        for(size_t j = 0; j < n_cols; j++){
            data[i][j] *= scalar;
        }
    }
    return *this;
}

Matrix& Matrix::iRowMulS(const Vector& other){

    assert(n_rows == other.size());

    for(size_t i = 0; i < n_rows; i++){
        double scalar = other.get(i);
        for(size_t j = 0; j < n_cols; j++){
            data[i][j] *= scalar;
        }
    }
    return *this;
}

Matrix& Matrix::iRowDivS(const Vector& other){

    assert(n_rows == other.size());

    for(size_t i = 0; i < n_rows; i++){
        double scalar = other.get(i);
        for(size_t j = 0; j < n_cols; j++){
            data[i][j] /= scalar;
        }
    }
    return *this;
}

Matrix Matrix::covariance() const {

    Matrix mat(n_cols, n_cols);

    for(size_t i = 0; i < n_cols; i++){
        for(size_t j = 0; j < n_cols; j++){
            double sum = 0.0f;
            for(size_t k = 0; k < n_rows; k++){
                sum += data[k][i] * data[k][j];
            }
            mat.set(i, j, sum);
        }
    }
    return mat;
}

Matrix Matrix::transpose() const {

    Matrix transpose(n_cols, n_rows);

    for(size_t i = 0; i < n_rows; i++){
        for(size_t j = 0; j < n_cols; j++){
            transpose.set(j, i, data[i][j]);
        }
    }
    return transpose;
}

// --- CSR ---

void
CSRMatrix::construct(std::vector<std::tuple<double, size_t, size_t>> raw_data,
                     size_t n_rows_, size_t n_cols_) {
    n_rows = n_rows_;
    n_cols = n_cols_;

    sortRawDataByRow(raw_data);

    indptr.push_back(0);

    for (size_t i = 0; i < raw_data.size(); ++i) {
        double val = std::get<0>(raw_data[i]);
        size_t row_idx = std::get<1>(raw_data[i]);
        size_t col_idx = std::get<2>(raw_data[i]);

        data.push_back(val);
        indices.push_back(col_idx);

        while (indptr.size() <= row_idx) {
            indptr.push_back(i);
        }

    }

    while (indptr.size() <= n_rows) {
        indptr.push_back(data.size());
    }
}

CSRMatrix::CSRMatrix(
    std::vector<std::tuple<double, size_t, size_t>> raw_data,
    size_t n_rows_, size_t n_cols_) {

    construct(raw_data, n_rows_, n_cols_);
}

CSRMatrix::CSRMatrix(const Matrix& dense) {
    vector<tuple<double, size_t, size_t>> raw_data;

    for (size_t i = 0; i < dense.nRows(); ++i) {
        for (size_t j = 0; j < dense.nCols(); ++j) {
            if (dense.get(i, j) != 0) {
                raw_data.push_back(make_tuple(dense.get(i, j), i, j));
            }
        }
    }
    construct(raw_data, dense.nRows(), dense.nCols());
}

Vector
CSRMatrix::getRow(size_t row_idx) const {
    Vector ret(n_cols);

    size_t start_idx = indptr[row_idx];
    size_t end_idx = indptr[row_idx+1];

    for (size_t k = start_idx; k < end_idx; ++k) {
        size_t i = indices[k];
        double mat_val = data[k];
        ret.set(i, mat_val);
    }
    return ret;
}

Vector
CSRMatrix::mul(const Vector& vec) const {
    assert(n_cols == vec.size());

    Vector ret(vec.size());

    for (size_t i = 0; i < n_rows; ++i) {
        size_t start_idx = indptr[i];
        size_t end_idx = indptr[i+1];

        double temp = 0.0;

        for (size_t k = start_idx; k < end_idx; ++k) {
            size_t j = indices[k];

            temp += data[k] * vec.get(j);
        }
        ret.set(i, temp);
    }

    return ret;
};

Vector
CSRMatrix::mulPartial(const Vector& vec, size_t n_cols_) const {
    assert(n_cols_ == vec.size());

    Vector ret(vec.size());

    for (size_t i = 0; i < n_rows; ++i) {
        size_t start_idx = indptr[i];
        size_t end_idx = indptr[i+1];

        for (size_t k = start_idx; k < end_idx; ++k) {
            size_t j = indices[k];

            if (j >= n_cols_)
                break;

            ret.set(i, ret.get(i) + data[k] * vec.get(j));
        }
    }

    return ret;
};

bool
CSRMatrix::posInData(size_t i, size_t j, size_t& ret) const {
    assert(0 <= i); assert(i < n_rows);
    assert(0 <= j); assert(j < n_cols);

    size_t start = indptr[j];
    size_t end = indptr[j+1];

    for (size_t k = start; k < end; ++k) {
        if (indices[k] == i) {
            ret = k;
            return true;
        }
    }
    return false;
}

/*
 * Definition of friend functions of Vector class
 */
void vecDot(Vector& dst, const Vector& src1, const Vector& src2){

    size_t size = src1.size();

    assert(size == src2.size());

    dst.resize(size);
    for(size_t i = 0; i < size; i++){
        dst.data[i] = src1.data[i]*src2.data[i];
    }
}

void vecScalarMul(Vector& dst, const Vector& src, double scalar){

    size_t size = src.size();

    dst.resize(size);
    for(size_t i = 0; i < size; i++){
        dst.data[i] = src.data[i]*scalar;
    }
}

void vecSub(Vector& dst, const Vector& src1, const Vector& src2){

    size_t size = src1.size();

    assert(size == src2.size());

    dst.resize(size);
    for(size_t i = 0; i < size; i++){
        dst.data[i] = src1.data[i]-src2.data[i];
    }
}

void vecAdd(Vector& dst, const Vector& src1, const Vector& src2){

    size_t size = src1.size();

    assert(size == src2.size());

    dst.resize(size);
    for(size_t i = 0; i < size; i++){
        dst.data[i] = src1.data[i]+src2.data[i];
    }
}

void matVecMul(Vector& dst, const Matrix& mat, const Vector& vec) {

    size_t col = mat.nCols();
    size_t row = mat.nRows();

    assert(col == vec.size());

    dst.resize(mat.nRows());

    for (size_t i = 0; i < row; i++) {

        double sum = .0f;

        for (size_t j = 0; j < col; j++) {
            sum += mat.data[i][j] * vec.data[j];
        }
        dst.data[i] = sum;
    }
}

void matVecMulPartialT(
    Vector& dst, const Matrix& mat,
    const Vector& vec, size_t n_rows_) {

    assert(n_rows_ == vec.size());

    for (size_t j = 0; j < mat.n_cols; ++j) {
        double temp = .0f;
        for (size_t i = 0; i < n_rows_; ++i) {
            temp += mat.data[i][j] * vec.get(i);
        }
        dst.set(j, temp);
    }
}

double l2norm(const Vector& vec){

    size_t size = vec.size();
    double res = .0f;

    for (size_t i = 0; i < size; i++) {
        double num = vec.data[i];
        num *= num;
        res += num;
    }

    return sqrt(res);
}

void copyCol(Matrix& dst, const Matrix& src, size_t to, size_t from){

    size_t dst_col = dst.nCols();
    size_t src_col = src.nCols();
    size_t row = dst.nRows();

    assert(to < dst_col);
    assert(from < src_col);
    assert(row == src.nRows());

    for(size_t i = 0; i < row; i++){
        dst.data[i][to] = src.data[i][from];
    }
}

void copyRow(Matrix& dst, const Matrix& src, size_t to, size_t from){

    size_t dst_row = dst.nRows();
    size_t src_row = src.nRows();
    size_t col = dst.nCols();

    assert(to < dst_row);
    assert(from < src_row);
    assert(col == src.nCols());

    for(size_t j = 0; j < col; j++){
        dst.data[to][j] = src.data[from][j];
    }
}

void matMulRowCoef(Vector& dst, const Matrix& src1, const Matrix& src2, size_t row){

    size_t src1_col = src1.nCols();
    size_t src1_row = src1.nRows();
    size_t src2_col = src2.nCols();
    size_t src2_row = src2.nRows();

    assert(row < src2_row);
    assert(src1_col == src2_col);

    dst.resize(src1_row);
    for(size_t i = 0; i < src1_row; i++){
        double sum = .0f;
        for(size_t j = 0; j < src1_col; j++){
            sum += src1.data[i][j] * src2.data[row][j];
        }
        dst.data[i] = sum;
    }
}

/*
 * Definition of friend functions of CSRMatrix class
 */
void spMatVecMul(Vector& dst, const CSRMatrix& mat, const Vector& vec) {
    assert(dst.size() == vec.size());
    assert(mat.n_cols == vec.size());

    for (size_t i = 0; i < mat.n_rows; ++i) {
        size_t start_idx = mat.indptr[i];
        size_t end_idx = mat.indptr[i+1];

        double temp = 0.0;

        for (size_t k = start_idx; k < end_idx; ++k) {
            size_t j = mat.indices[k];

            temp += mat.data[k] * vec.get(j);
        }
        dst.set(i, temp);
    }
}


/*
 * Definition of helper functions used by CSRMatrix
 */
void sortRawDataByRow(std::vector<std::tuple<double, size_t, size_t>>& raw_data) {
    std::sort(
        raw_data.begin(), raw_data.end(),
        [](const std::tuple<double, size_t, size_t> &left,
           const std::tuple<double, size_t, size_t> &right) {
            if (std::get<1>(left) < std::get<1>(right)) {
                return true;
            } else if (std::get<1>(left) > std::get<1>(right)) {
                return false;
            } else {
                return std::get<2>(left) < std::get<2>(right);
            }
        });
}
