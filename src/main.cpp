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
gmres(const Matrix& A,
      const Vector& b,
      size_t m, double tol, size_t maxit) {

    char buf[1024];

    size_t dim = A.nCols();
    size_t nit = 0;
    size_t innit = 0;
    size_t outnit = 0;

    Matrix H = Matrix(m+1, m);
    Matrix Z = Matrix(m, dim);
    Matrix V = Matrix(m+1, dim);
    Vector x(dim);
    Vector x0(dim);

    double tKrylov = .0f;
    double tLLS = .0f;
    double tStart;

    assert(dim == b.size());

    // Use trivial preconditioner for now
    // auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    while (nit < maxit) {

        Vector r0 = b.sub(A.mul(x0));
        double beta = r0.norm2();
        V.setRow(0, r0.mulS(1.0 / beta));

        innit = 0;
        
        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {

            tStart = CycleTimer::currentSeconds();

            // Z.setCol(j, Prcnd.mul(V.getCol(j)));
            Z.setRow(j, V.getRow(j));

            Vector w = A.mul(Z.getRow(j));

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

            //x = x0.add(Z.mulPartial(y, j+1));
            x = x0.add(Z.mulPartialT(y, j+1));

            double res_norm = A.mul(x).sub(b).norm2();

            nit++;
            innit++;
            
            tLLS += CycleTimer::currentSeconds() - tStart;

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
sparseGmres(const CSRMatrix& A,
         const Vector& b,
         size_t m, double tol, size_t maxit) {

    char buf[1024];

    size_t dim = A.nCols();
    size_t nit = 0;
    size_t innit = 0;
    size_t outnit = 0;

    Matrix H = Matrix(m+1, m);
    Matrix Z = Matrix(m, dim);
    Matrix V = Matrix(m+1, dim);
    Vector x0(dim);
    Vector x(dim);

    double tKrylov = .0f;
    double tLLS = .0f;
    double tStart;

    assert(dim == b.size());

    // Use trivial preconditioner for now
    //auto Prcnd = identitySparseMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;


    while (nit < maxit) {


        Vector r0 = b.sub(A.mul(x0));
        double beta = r0.norm2();
        V.setRow(0, r0.mulS(1.0 / beta));

        innit = 0;
        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {

            tStart = CycleTimer::currentSeconds();

            //Z.setRow(j, Prcnd.mul(V.getRow(j)));
            Z.setRow(j, V.getRow(j));

            Vector w = A.mul(Z.getRow(j));

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

            x = x0.add(Z.mulPartialT(y, j+1));

            double res_norm = A.mul(x).sub(b).norm2();

            nit++;
            innit++;

            tLLS += CycleTimer::currentSeconds() - tStart;

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
ompGmres(const Matrix& A,
         const Vector& b,
         size_t m, double tol, size_t maxit) {

    char buf[1024];

    size_t dim = A.nCols();
    size_t nit = 0;
    size_t innit = 0;
    size_t outnit = 0;

    Matrix H = Matrix(m+1, m);
    Matrix Z = Matrix(dim, m);
    Matrix V = Matrix(dim, m+1);
    Vector x(dim);
    Vector x0(dim);
    Vector r0;
    Vector temp;

    double tKrylov = .0f;
    double tLLS = .0f;
    double tStart;

    assert(dim == b.size());

    // Use trivial preconditioner for now
    // auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    while (nit < maxit) {

        matVecMul(temp, A, x0);
        vecSub(r0, b, temp);
        double beta = r0.norm2();
        vecScalarMul(temp, r0, 1.0/beta);
        V.setRow(0, temp);

        innit = 0;
        
        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {

            tStart = CycleTimer::currentSeconds();

            // Z.setCol(j, Prcnd.mul(V.getCol(j)));
            Z.setCol(j, V.getCol(j));

            Vector w = A.mul(Z.getCol(j));

            for (size_t i = 0; i < j; i++) {
                Vector v = V.getCol(i);
                H.set(i, j, w.dotV(v));
                w.isub(v.mulS(H.get(i, j)));
            }

            H.set(j+1, j, w.norm2());
            V.setCol(j+1, w.mulS(1.0 / H.get(j+1, j)));


            tKrylov += CycleTimer::currentSeconds() - tStart;
            tStart = CycleTimer::currentSeconds();

            Vector y = leastSquareWithQR(H, j+1, beta);

            x = x0.add(Z.mulPartial(y, j+1));

            double res_norm = A.mul(x).sub(b).norm2();

            nit++;
            innit++;
            
            tLLS += CycleTimer::currentSeconds() - tStart;

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

    // experiment on dense matrix representation
    start_time = CycleTimer::currentSeconds();
    gmres(A, b, m, tol, maxit);
    end_time = CycleTimer::currentSeconds();
    
    sprintf(buf, "[%.3f] ms in total (Dense)\n\n", 
            (end_time - start_time) * 1000);
    cout << buf;

    // experiment on scr sparse matrix representation
    start_time = CycleTimer::currentSeconds();
    sparseGmres(A_csr, b, m, tol, maxit);
    end_time = CycleTimer::currentSeconds();
    
    sprintf(buf, "[%.3f] ms in total (CSR Sparse GMRES) \n\n", 
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
