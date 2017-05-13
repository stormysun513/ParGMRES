#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>
#include <cmath>
#include <string>

#include "CycleTimer.h"

#include "utils.h"

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100000

using namespace std;



Vector
sparseGmres(
    const CSRMatrix& A, const Vector& b,
    size_t m, double tol, size_t maxit) {

    char buf[1024];

    size_t dim = A.nCols();
    size_t nit = 0;
    size_t innit = 0;
    size_t outnit = 0;

    Matrix H = Matrix(m+1, m);
    Matrix V = Matrix(m+1, dim);
    Vector x(dim);
    Vector x0(dim);
    Vector r0;
    Vector temp(dim);
    Vector w(dim);

    double tKrylov = .0f;
    double tLLS = .0f;
    double tGetRes = .0f;
    double tStart;

    assert(dim == b.size());

    // Use trivial preconditioner for now
    // auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    while (nit < maxit) {

        spMatVecMul(temp, A, x0);
        vecSub(r0, b, temp);
        double beta = r0.norm2();
        vecScalarMul(temp, r0, 1.0/beta);
        V.setRow(0, temp);

        innit = 0;

        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {

            tStart = CycleTimer::currentSeconds();

            spMatVecMul(w, A, V.getRow(j));

            for (size_t i = 0; i < j; i++) {
                Vector v = V.getRow(i);
                H.set(i, j, w.dotV(v));
                w.isub(v.mulS(H.get(i, j)));
            }

            H.set(j+1, j, w.norm2());
            V.setRow(j+1, w.mulS(1.0 / H.get(j+1, j)));


            tKrylov += CycleTimer::currentSeconds() - tStart;
            tStart = CycleTimer::currentSeconds();

            Vector y = leastSquareWithQR(H, j+1, beta);

            tLLS += CycleTimer::currentSeconds() - tStart;
            tStart = CycleTimer::currentSeconds();

            matVecMulPartialT(temp, V, y, j+1);
            x = x0.add(temp);

            spMatVecMul(temp, A, x);
            double res_norm = temp.sub(b).norm2();

            tGetRes += CycleTimer::currentSeconds() - tStart;

            nit++;
            innit++;

            if (res_norm < tol * b.norm2()) {
                cout << "FGMRES converged to relative tolerance: "
                     << res_norm / b.norm2()
                     << " at iteration "
                     << nit
                     << " (out: "
                     << outnit
                     << ", in: "
                     << innit
                     << ")"
                     << endl;

                sprintf(buf, "[%.3f] ms in Krylov \n", tKrylov * 1000);
                cout << buf;
                sprintf(buf, "[%.3f] ms in LLS \n", tLLS * 1000);
                cout << buf;
                sprintf(buf, "[%.3f] ms in Get Residual \n", tGetRes * 1000);
                cout << buf;

                return x;
            }
        }
        x0 = x;
        outnit++;
    }

    spMatVecMul(temp, A, x0);
    double res_norm = temp.sub(b).norm2();
    cout << "FGMRES is not converged: "
         << res_norm / b.norm2()
         << endl;
    return x0;
}

void runExp(const string& mat_name) {

    int m = 100;                    // m may be adaptive to input size
    int maxit = 100000;
    double tol = 1e-6;
    double start_time;
    double end_time;
    char buf[1024];

    Matrix A = loadMTXToMatrix(mat_name);
    assert(A.nCols() == A.nRows());

    CSRMatrix A_csr(A);
    Vector b = randUnitUniformVector(A.nCols());

    cout << "A: " << mat_name << " "
         << A.nRows() << "x" << A.nCols() << endl;
    cout << "m=" << m << ", tol=" << tol << ", maxit=" << maxit << endl;
    cout << endl;

    // experiment on CSR + OMP
    start_time = CycleTimer::currentSeconds();
    sparseGmres(A_csr, b, m, tol, maxit);
    end_time = CycleTimer::currentSeconds();

    sprintf(buf, "[%.3f] ms in total (Sparse w/o OMP) \n\n",
            (end_time - start_time) * 1000);
    cout << buf;
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " matrix_filename" << endl;
        return 1;
    }
    runExp(argv[1]);

    return 0;
}
