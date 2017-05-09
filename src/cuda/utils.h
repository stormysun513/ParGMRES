#ifndef __UTILS_H__
#define __UTILS_H__

#include "mtxvec.h"

/* DEBUG helper function */
void printVector(const Vector& vec);
void printMatrix(const Matrix& mat, size_t row_start, size_t row_end,
                 size_t col_start, size_t col_end);

/* vector and matrix generation function */
Matrix identityMatrix(size_t size);

/* QR */
Vector leastSquareWithQR(const Matrix& H, size_t size, float beta);
void householderQR(const Matrix& H, size_t size, Matrix& Q, Matrix& R);
void qrSolve(const Matrix& Q, const Matrix& R, const Vector& b, Vector& sol);

#endif  // UTILS_H
