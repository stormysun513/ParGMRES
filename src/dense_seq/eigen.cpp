#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <algorithm>

#include "../Eigen/Dense"

#include "utils.h"

#define MAX_SVD_ITERATION   500
#define SVD_EPSILON         1e-3

static Vector svdOneDimPower(const Matrix& A){

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

static void svdPowerMethod(Matrix& U, Matrix& V, Vector& sigma, const Matrix& A){
    
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

        Vector v = svdOneDimPower(tmp);
        Vector u_unnormalized = A.mul(v);
        double lambda = u_unnormalized.norm2();

        sigma.set(count, lambda);
        V.setRow(count, v);
        U.setRow(count, u_unnormalized.normalize());

        count++;
    }
}

void svdJacobiMethod(Matrix& U, Matrix& V, Vector& sigma, const Matrix& A){
    
    size_t n = A.nCols();
    size_t iters = 0;

    U = A;
    V = identityMatrix(n);
    sigma.resize(n);

    while(iters < MAX_SVD_ITERATION){
 
        double converge = 0;

        // O(n^3)
        for(size_t j = 1; j < n; j++){
            for(size_t k = 0; k <= j-1; k++){
            
                Vector ck = U.getCol(k);
                Vector cj = U.getCol(j);

                double alpha = ck.dotV(ck);
                double beta = cj.dotV(cj); 
                double gamma = ck.dotV(cj);
                double tan;

                converge = std::fmax(converge, (std::fabs(gamma)/std::sqrt(alpha*beta)));


                if(gamma != 0){
                    double zeta = (beta-alpha)/(2*gamma);
                    tan = 1.0/(std::fabs(zeta)+std::sqrt(1 + zeta*zeta));
                    tan = (zeta > 0) ? tan : -1*tan;
                }
                else{
                    tan = 0;
                }
                
                double cos = 1.0/std::sqrt(tan*tan+1);
                double sin = cos * tan;

                //Vector tk = U.getCol(k);
                U.setCol(k, ck.mulS(cos).sub(U.getCol(j).mulS(sin)));
                U.setCol(j, ck.mulS(sin).add(U.getCol(j).mulS(cos)));

                ck = V.getCol(k);
                V.setCol(k, ck.mulS(cos).sub(V.getCol(j).mulS(sin)));
                V.setCol(j, ck.mulS(sin).add(V.getCol(j).mulS(cos)));
            }
        }

        if(converge < SVD_EPSILON){
            break;
        }

        iters++;
    }

    for(size_t i = 0; i < n; i++){
        sigma.set(i, U.getCol(i).norm2());
    }
    U = U.transpose().iRowDivS(sigma);
}

Vector leastSquareWithJacobi(const Matrix& H, size_t size, double beta){

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

    svdJacobiMethod(U, V, sigma, A);

    size_t length = sigma.size();
    for(size_t i = 0; i < length; i++){
        sigma.set(i, 1.0/sigma.get(i));
    }

    Matrix UT = U.iRowMulS(sigma);
    Vector y = V.mul(U.mul(b));

    return y;
}

Vector leastSquareWithPowerMethod(const Matrix& H, size_t size, double beta){

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

    svdPowerMethod(U, V, sigma, A);


    size_t length = sigma.size();
    for(size_t i = 0; i < length; i++){
        sigma.set(i, 1.0/sigma.get(i));
    }

    Matrix UT = U.iRowMulS(sigma);
    Matrix UV = V.transpose();

    Vector y = UV.mul(UT.mul(b));

    return y;
}


Vector leastSquareWithEigen(const Matrix& H, size_t size, double beta) {

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
