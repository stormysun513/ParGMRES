#include "utils.hpp"
#include "loadmtx.hpp"
#include "mtxvec.hpp"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "Eigen/Dense"

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100

using namespace std;

Matrix
identityMatrix(size_t size) {
    Matrix mat(size, size);

    for (size_t i = 0; i < size; i++) {
        mat.set(i, i, 1.0);
    }
    return mat;
}

void
printMatrix(Matrix& mat,
            size_t row_start, size_t row_end,
            size_t col_start, size_t col_end) {
    size_t m = mat.nRows();
    size_t n = mat.nCols();

    assert(row_start >= 0);
    assert(col_start >= 0);
    assert(row_end <= m);
    assert(col_end <= n);

    char buf[1024];
    sprintf(buf, "");

    for (size_t i = row_start; i < row_end; i++) {
        for (size_t j = col_start; j < col_end; j++) {
            if (j == col_start)
                sprintf(buf, "%s[", buf);

            sprintf(buf, "%s%.4f", buf, mat.get(i, j));

            if (j == col_end-1)
                sprintf(buf, "%s]", buf);
            else
                sprintf(buf, "%s, ", buf);
        }
        sprintf(buf, "%s\n", buf);
    }

    cout << buf;
}

Vector
leastSquare(Matrix& H, size_t size, double beta) {
    Eigen::MatrixXf A(size+1, size);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(size+1);
    Eigen::VectorXf y_(size);
    Vector y(size);

    b(0) = beta;
    for (int i = 0; i < size+1; i++) {
        for (int j = 0; j < size; j++) {
            A(i, j) = H.get(i, j);
        }
    }

    y_ = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    for (int i = 0; i < size; i++) {
        y.set(i, y_(i));
    }

    return y;
}

Vector
gmres(const Matrix& A,
      const Vector& b,
      size_t m, double tol, size_t maxit) {

    size_t dim = A.nCols();
    size_t nit = 0;
    Vector x0(dim);

    assert(dim == b.size());

    // Use trivial preconditioner for now
    auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    Vector x(dim);
    while (nit < maxit) {
        Matrix H = Matrix(m+1, m);
        Matrix Z = Matrix(dim, m);
        Matrix V = Matrix(dim, m+1);

        Vector r0 = b.sub(A.mul(x0));
        double beta = r0.norm2();
        V.setCol(0, r0.mulS(1.0 / beta));

        // Generate krylov subspace
        for(int j = 0; j < m; j++) {
            // Z[:, j] = P * V[:, j]
            Z.setCol(j, Prcnd.mul(V.getCol(j)));

            // w = A * Z[:, j]
            Vector w = A.mul(Z.getCol(j));

            for (size_t i = 0; i < j; i++) {
                Vector v = V.getCol(i);
                H.set(i, j, w.dotV(v));
                w.isub(v.mulS(H.get(i, j)));
            }

            H.set(j+1, j, w.norm2());
            V.setCol(j+1, w.mulS(1.0 / H.get(j+1, j)));

            // printMatrix(H, 0, j+1, 0, j);
            // cout << "---" << endl;

            Vector y = leastSquare(H, j+1, beta);
            x = x0.add(Z.mulPartial(y, j+1));

            double res_norm = A.mul(x).sub(b).norm2();

            if (res_norm < tol * b.norm2()) {
                cout << "FGMRES converged to relative tolerance: "
                     << res_norm / b.norm2()
                     << " at iteration "
                     << nit << endl;
                return x;
            }
        }

        x0 = x;
        nit++;
    }
    return x;
}

int main(int argc, char *argv[])
{
    int m = 100;
    int maxit = 100;
    double tol = 1e-3;

    Matrix A = loadMTXFile("../data/cage4.mtx");
    Vector b = Vector(A.nCols());

    cout << "A = ../data/cage4.mtx" << endl;
    b.set(0, 1.0);

    cout << "m=" << m << ", tol=" << tol << ", maxit=" << maxit << endl;
    gmres(A, b, m, tol, maxit);

    tol = 1e-9;
    cout << "m=" << m << ", tol=" << tol << ", maxit=" << maxit << endl;
    gmres(A, b, m, tol, maxit);

    return 0;
}
