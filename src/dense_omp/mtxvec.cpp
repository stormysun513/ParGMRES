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

void matVecMul(Vector& dst, const Matrix& mat, const Vector& vec){

    size_t col = mat.nCols();
    size_t row = mat.nRows();

    assert(col == vec.size());

    dst.resize(mat.nRows());

    if (col < OMP_NN_BOUND) {
        for (size_t i = 0; i < row; i++) {

            double sum = .0f;

            for (size_t j = 0; j < col; j++) {
                sum += mat.data[i][j] * vec.data[j];
            }
            dst.data[i] = sum;
        }
    } else {
#pragma omp parallel for schedule(static, 64)
        for (size_t i = 0; i < row; i++) {

            double sum = .0f;

            for (size_t j = 0; j < col; j++) {
                sum += mat.data[i][j] * vec.data[j];
            }
            dst.data[i] = sum;
        }
    }
}

void matVecMulPartialT(
    Vector& dst, const Matrix& mat,
    const Vector& vec, size_t n_rows_) {

    assert(n_rows_ == vec.size());

    if (mat.n_cols < OMP_NN_BOUND) {
        for (size_t j = 0; j < mat.n_cols; ++j) {
            double temp = .0f;
            for (size_t i = 0; i < n_rows_; ++i) {
                temp += mat.data[i][j] * vec.get(i);
            }
            dst.set(j, temp);
        }
    } else {
#pragma omp parallel for schedule(static, 64)
        for (size_t j = 0; j < mat.n_cols; ++j) {
            double temp = .0f;
            for (size_t i = 0; i < n_rows_; ++i) {
                temp += mat.data[i][j] * vec.get(i);
            }
            dst.set(j, temp);
        }
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
