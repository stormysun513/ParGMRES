#ifndef __UTILS_H__
#define __UTILS_H__

#include "mtxvec.h"

/* DEBUG helper function */
void printVector(const Vector& vec);
void printMatrix(Matrix& mat, size_t row_start, size_t row_end,
            size_t col_start, size_t col_end);

/* vector and matrix generation function */
Vector randUniformVector(size_t n);
Vector randUnitUniformVector(size_t n);
SparseMatrix identitySparseMatrix(size_t size);
Matrix identityMatrix(size_t size);

/* algorithms */
Vector leastSquare(const Matrix& H, size_t size, double beta);
Vector leastSquareWithBeta(const Matrix& H, size_t size, double beta);
void svd(Matrix& U, Matrix& V, Vector& sigma, const Matrix& A);

/* mtx file IO */
Vector loadMTXToVector(const std::string&);
Matrix loadMTXToMatrix(const std::string&);
void writeVecToMTXFile(const std::string&, const Vector&);
void writeMatToMTXFile(const std::string&, const Matrix&);

#endif  // UTILS_H
