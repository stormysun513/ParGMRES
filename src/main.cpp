#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>
#include <cmath>
#include <string>
#include <random>

#include "CycleTimer.h"

#include "utils.h"

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100000

using namespace std;


Vector
gmres(const Matrix& A,
      const Vector& b,
      size_t m, double tol, size_t maxit) {

    size_t dim = A.nCols();
    size_t nit = 0;
    size_t innit = 0;
    size_t outnit = 0;
    Vector x0(dim);

    assert(dim == b.size());

    // Use trivial preconditioner for now
    auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    while (nit < maxit) {
        Matrix H = Matrix(m+1, m);
        Matrix Z = Matrix(dim, m);
        Matrix V = Matrix(dim, m+1);
        Vector x(dim);

        Vector r0 = b.sub(A.mul(x0));
        double beta = r0.norm2();
        V.setCol(0, r0.mulS(1.0 / beta));

        innit = 0;
        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {
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

            Vector y = leastSquare(H, j+1, beta);
            //Vector y = leastSquareWithBeta(H, j+1, beta);
            //
            x = x0.add(Z.mulPartial(y, j+1));

            double res_norm = A.mul(x).sub(b).norm2();

            nit++;
            innit++;
            if (res_norm < tol * b.norm2()) {
                cout << "FGMRES converged to relative tolerance: "
                     << res_norm / b.norm2()
                     << " at iteration "
                     << nit
                     << "(out: "
                     << outnit
                     << ", in: "
                     << innit
                     << ")"
                     << endl;
                return x;
            }
        }
        x0 = x;
        outnit++;
    }

    double res_norm = A.mul(x0).sub(b).norm2();
    cout << "FGMRES is not converged: "
         << res_norm / b.norm2()
         << endl;
    return x0;
}

Vector
sparseGmres(const SparseMatrix& A,
            const Vector& b,
            size_t m, double tol, size_t maxit) {

    size_t dim = A.nCols();
    size_t nit = 0;
    size_t innit = 0;
    size_t outnit = 0;
    Vector x0(dim);

    assert(dim == b.size());

    // Use trivial preconditioner for now
    auto Prcnd = identitySparseMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    while (nit < maxit) {
        Matrix H = Matrix(m+1, m);

        Matrix Z = Matrix(m, dim);
        Matrix V = Matrix(m+1, dim);
        Vector x(dim);

        Vector r0 = b.sub(A.mul(x0));
        double beta = r0.norm2();
        V.setRow(0, r0.mulS(1.0 / beta));

        innit = 0;
        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {
            // Z[j, :] = P * V[j, :]
            Z.setRow(j, Prcnd.mul(V.getRow(j)));

            // w = A * Z[j, :]
            Vector w = A.mul(Z.getRow(j));

            for (size_t i = 0; i < j; i++) {
                Vector v = V.getRow(i);
                H.set(i, j, w.dotV(v));
                w.isub(v.mulS(H.get(i, j)));
            }

            H.set(j+1, j, w.norm2());
            V.setRow(j+1, w.mulS(1.0 / H.get(j+1, j)));

            Vector y = leastSquare(H, j+1, beta);
            //Vector y = leastSquareWithBeta(H, j+1, beta);

            x = x0.add(Z.mulPartialT(y, j+1));

            double res_norm = A.mul(x).sub(b).norm2();

            nit++;
            innit++;
            if (res_norm < tol * b.norm2()) {
                cout << "FGMRES converged to relative tolerance: "
                     << res_norm / b.norm2()
                     << " at iteration "
                     << nit
                     << "(out: "
                     << outnit
                     << ", in: "
                     << innit
                     << ")"
                     << endl;
                return x;
            }
        }
        x0 = x;
        outnit++;
    }

    double res_norm = A.mul(x0).sub(b).norm2();
    cout << "FGMRES is not converged: "
         << res_norm / b.norm2()
         << endl;
    return x0;
}

void runExp(const string& mat_name) {
    int m = 100;
    int maxit = 1000000;
    double tol = 1e-6;
    double start_time;
    double end_time;
    char buf[1024];

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    Matrix A = loadMTXToMatrix(mat_name);
    SparseMatrix A_csc(A);
    Vector b = Vector(A.nCols());

    for (size_t i = 0; i < A.nCols(); ++i) {
        b.set(i, distribution(generator));
    }

    cout << "A: " << mat_name << " "
         << A.nRows() << "x" << A.nCols() << endl;

    cout << "m=" << m << ", tol=" << tol << ", maxit=" << maxit << endl;
    start_time = CycleTimer::currentSeconds();
    gmres(A, b, m, tol, maxit);
    end_time = CycleTimer::currentSeconds();
    sprintf(buf, "[%.3f] ms\n\n", (end_time - start_time) * 1000);
    cout << buf;

    start_time = CycleTimer::currentSeconds();
    sparseGmres(A_csc, b, m, tol, maxit);
    end_time = CycleTimer::currentSeconds();
    sprintf(buf, "[%.3f] ms (Sparse GMRES) \n\n", (end_time - start_time) * 1000);
    cout << buf;
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " matrix_filename" << endl;
        return 0;
    }
    runExp(argv[1]);

    return 0;
}
