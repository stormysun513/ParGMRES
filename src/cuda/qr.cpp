#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <algorithm>
#include <cassert>

#include "utils.h"

static float partialNorm(const Matrix& R, size_t j);
static Matrix getwwTR(const Matrix& R, float* w, size_t j, float tau);
static Matrix getQwwT(const Matrix& Q, float* w, size_t j, float tau);

Vector leastSquareWithQR(const Matrix& H, size_t size, float beta) {
    Matrix Q;
    Matrix R;
    Vector sol(size);
    Vector b(size+1);
    b.set(0, beta);

    householderQR(H, size, Q, R);
    qrSolve(Q, R, b, sol);

    return sol;
}

// Perform QR on H[:(size+1), :size]
// usage:
//     std::default_random_engine generator;
//     std::uniform_real_distribution<float> distribution(-1.0, 1.0);

//     Matrix A(5, 4);
//     A.set(0, 0, distribution(generator));
//     A.set(1, 0, distribution(generator));
//     A.set(0, 1, distribution(generator));
//     A.set(1, 1, distribution(generator));
//     A.set(2, 1, distribution(generator));
//     A.set(0, 2, distribution(generator));
//     A.set(1, 2, distribution(generator));
//     A.set(2, 2, distribution(generator));
//     A.set(3, 2, distribution(generator));
//     A.set(0, 3, distribution(generator));
//     A.set(1, 3, distribution(generator));
//     A.set(2, 3, distribution(generator));
//     A.set(3, 3, distribution(generator));
//     A.set(4, 3, distribution(generator));

//     printMatrix(A, 0, 5, 0, 4);
//     Matrix Q;
//     Matrix R;
//     householderQR(A, 4, Q, R);
//     printMatrix(Q, 0, 5, 0, 5);
//     printMatrix(R, 0, 5, 0, 4);
void householderQR(const Matrix& H, size_t size, Matrix& Q, Matrix& R) {
    size_t m = size+1;
    size_t n = size;
    Q = identityMatrix(m);
    R = Matrix(m, n);

    // copy R = H
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            R.set(i, j, H.get(i, j));
        }
    }

    for (size_t j = 0; j < n; ++j) {
        // nomrx = norm(R[j:, j])
        float normx = partialNorm(R, j);
        float s = (R.get(j, j) > 0 ? -1.0 : 1.0);
        float u1 = R.get(j, j) - s * normx;

        float w[m-j];
        w[0] = 1.0;
        for(size_t i = j+1; i < m; ++i) {
            w[i-j] = R.get(i, j) / u1;
        }

        float tau = -1 * s * u1 / normx;

        // R(j:end,:) = R(j:end,:)-(tau*w)*(w’*R(j:end,:));
        Matrix sndTerm = getwwTR(R, w, j, tau);
        for (size_t ii = j; ii < m; ++ii) {
            for (size_t jj = 0; jj < n; ++jj) {
                R.set(ii, jj, R.get(ii, jj)- sndTerm.get(ii - j, jj));
            }
        }

        // Q(:,j:end) = Q(:,j:end)-(Q(:,j:end)*w)*(tau*w)’;
        sndTerm = getQwwT(Q, w, j, tau);
        for (size_t ii = 0; ii < m; ++ii) {
            for (size_t jj = j; jj < n; ++jj) {
                Q.set(ii, jj, Q.get(ii, jj)- sndTerm.get(ii, jj - j));
            }
        }
    }
}

// Usage:
//     Vector b(5);
//     b.set(0, 10.0);
//     b.set(1, 0.0);
//     b.set(2, -3.0);
//     b.set(3, 4.0);
//     b.set(4, -1.0);
//     Vector sol(4);
//     qrSolve(Q, R, b, sol);
void qrSolve(const Matrix& Q, const Matrix& R, const Vector& b, Vector& sol) {
    // Solve Ax = b ==> Solve Rx = Q^T b

    assert(Q.nCols() == Q.nRows());
    assert(Q.nRows() == b.size());
    assert(sol.size() == R.nCols());

    // Q^T b, shape == b.shape
    Vector QTb(b.size());

    for (size_t i = 0; i < Q.nCols(); ++i) {
        float temp = 0.0;

        for (size_t j = 0; j < Q.nRows(); ++j) {
            temp += Q.get(j, i) * b.get(j);
        }

        QTb.set(i, temp);
    }

    // Back-substitute

    for (size_t i = 0; i < R.nCols(); ++i) {
        size_t row_idx = R.nCols() - 1 - i;
        float res = QTb.get(row_idx);

        for (size_t j = row_idx+1; j < R.nCols(); ++j) {
            res -= sol.get(j) * R.get(row_idx, j);
        }

        sol.set(row_idx,  res / R.get(row_idx, row_idx));
    }
}

static float partialNorm(const Matrix& R, size_t j) {
    float ret = 0.0;
    for (size_t i = j; i < R.nRows(); ++i) {
        ret += R.get(i, j) * R.get(i, j);
    }

    return sqrt(ret);
}


// w = len x 1,
// w^T = 1 x len
// R = len x (n_cols)
// w^T dot R = 1 x n_cols

static Matrix getwwTR(const Matrix& R, float* w, size_t j, float tau) {
    size_t len = (R.nRows() - j);
    // (tau*w)*(w’*R(j:end,:))
    float wTR[R.nCols()];
    Matrix wwTR(len, R.nCols());

    // wTR = w^T dot R, shape: 1 x n_cols
    for (size_t k = 0; k < R.nCols(); ++k) {
        float temp = 0.0;

        for (size_t i = 0; i < len; ++i) {
            temp += R.get(j+i, k) * w[i];
        }

        wTR[k] = temp;
    }

    // w dot wTR, shape: len x n_cols
    for (size_t i = 0; i < len; ++i) {
        for (size_t j = 0; j < R.nCols(); ++j) {
            wwTR.set(i, j, w[i] * wTR[j] * tau);
        }
    }

    return wwTR;
}

// w = len x 1,
// w^T = 1 x len
// Q = n_rows x len
// Q w = n_rows x 1
// Q w w^T = n_rows x len

static Matrix getQwwT(const Matrix& Q, float* w, size_t j, float tau) {
    size_t len = (Q.nRows() - j);
    float Qw[Q.nRows()];
    Matrix QwwT(Q.nRows(), len);

    // Qw = Q x w, shape: n_rows x 1
    for (size_t k = 0; k < Q.nRows(); ++k) {
        float temp = 0.0;

        for (size_t i = 0; i < len; ++i) {
            temp += Q.get(k, j+i) * w[i];
        }

        Qw[k] = temp;
    }

    // Q w w^T tau, shape: n_rows x len
    for (size_t i = 0; i < Q.nRows(); ++i) {
        for (size_t j = 0; j < len; ++j) {
            QwwT.set(i, j, Qw[i] * w[j] * tau);
        }
    }

    return QwwT;
}
