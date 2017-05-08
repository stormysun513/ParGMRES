#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <algorithm>

#include "utils.h"

#define MAX_SVD_ITERATION   500
#define SVD_EPSILON         1e-3

Vector leastSquareWithQR(const Matrix& H, size_t size, double beta) {
    Matrix Q;
    Matrix R;
    Vector sol(size);
    Vector b(size+1);
    b.set(0, beta);

    householderQR(H, size, Q, R);
    qrSolve(Q, R, b, sol);

    return sol;
}
