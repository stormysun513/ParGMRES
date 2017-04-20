#include <iostream>
#include <cmath>
#include <random>

#include "Eigen/Dense"

#include "utils.h"

#define MAX_SVD_ITERATION   1000
#define SVD_EPSILON         1e-3

static Vector svdOneDim(const Matrix& A){

    size_t n = A.nCols();

    Matrix B = A.covariance();
    Vector last(n);
    Vector curr = randUnitUniformVector(n);

    size_t iters = 0;
    while(iters <= MAX_SVD_ITERATION){
        last.copy(curr);
        curr = B.mul(curr).normalize();

        double cos = curr.dotV(last);
        if(cos > 1.0 - SVD_EPSILON){
            return curr;
        }
        iters++;
    }
    return curr;
}

void svd(Matrix& U, Matrix& V, Vector& sigma, const Matrix& A){

    size_t m = A.nRows();
    size_t n = A.nCols();
    size_t k = std::min(m, n);

    Matrix tmp = A;
    size_t count = 0;

    U.resize(k, m);
    V.resize(n, n);
    sigma.resize(k);
    for(size_t i = 0; i < k; i++){

        size_t j = count - 1;

        if(j < count){
            Vector v = V.getRow(j);
            Vector u = U.getRow(j);
            tmp.isub(u.crossV(v).imulS(sigma.get(j)));
        }

        Vector v = svdOneDim(tmp);
        Vector u_unnormalized = A.mul(v);
        double lambda = u_unnormalized.norm2();

        sigma.set(count, lambda);
        V.setRow(count, v);
        U.setRow(count, u_unnormalized.normalize());

        count++;
    }
}

Vector leastSquareWithBeta(const Matrix& H, size_t size, double beta){

    Matrix U, V;
    Matrix A(size+1, size);
    Vector sigma;
    Vector b(size+1);

    b.set(0, beta);
    for(size_t i = 1; i < size+1; i++){
        b.set(i, 0);
    }
    for(size_t i = 0; i < size+1; i++){
        for(size_t j = 0; j < size; j++){
            A.set(i, j, H.get(i, j));
        }
    }

    svd(U, V, sigma, A);


    size_t length = sigma.size();
    for(size_t i = 0; i < length; i++){
        sigma.set(i, 1.0/sigma.get(i));
    }

    Matrix UT = U.iRowMulS(sigma);
    Matrix UV = V.transpose();

    Vector y = UV.mul(UT.mul(b));

    return y;
}


Vector leastSquare(const Matrix& H, size_t size, double beta) {
    Eigen::MatrixXf A(size+1, size);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(size+1);
    Eigen::VectorXf y_(size);
    Vector y(size);

    b(0) = beta;
    for (size_t i = 0; i < size+1; i++) {
        for (size_t j = 0; j < size; j++) {
            A(i, j) = H.get(i, j);
        }
    }

    y_ = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    for (size_t i = 0; i < size; i++) {
        y.set(i, y_(i));
    }
    return y;
}
