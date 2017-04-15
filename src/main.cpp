#include "main.hpp"
#include "utils.hpp"
#include "loadmtx.hpp"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "Eigen/Dense"

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100

using namespace std;

Vector
mvAx(const vector<vector<double>>& A, const Vector& x) {
    Vector b(A.size());

    for (size_t i = 0; i < A.size(); ++i) {
        b.set(i, 0);
        for (size_t j = 0; j < A[0].size(); ++j) {
            b.set(i, b.get(i) + A[i][j] * x.get(j));
        }
    }
    return b;
}

// Like mvAx, but special case for Z(:, 1:j) * y
Vector
mvZy(const vector<vector<double>>& Z, const Vector& y, size_t j) {
    assert(j == y.size());
    Vector ret(Z.size());

    for (size_t i = 0; i < Z.size(); ++i) {
        ret.set(i, 0);
        for (size_t k = 0; k < j; ++k) {
            ret.set(i, ret.get(i) + Z[i][k] * y.get(k));
        }
    }
    return ret;
}

vector<vector<double>>
generateMatrix(int m, int n) {
    vector<vector<double>> matrix(m);

    for (size_t i = 0; i < m; i++) {
        matrix[i] = vector<double>(n);
    }
    return matrix;
}

vector<vector<double>>
identityMatrix(size_t size) {
    vector<vector<double>> matrix(size);

    for (size_t i = 0; i < size; i++) {
        matrix[i] = vector<double>(size);
        matrix[i][i] = 1.0;
    }

    return matrix;
}

void
setRow(vector<vector<double>>& mat, const Vector& vec, int row) {

    size_t length = vec.size();

    assert(mat.size() > row);
    assert(mat[0].size() == length);

    for (int i = 0; i < length; i++) {
        mat[row][i] = vec.get(i);
    }
}

void
setCol(vector<vector<double>>& mat, const Vector& vec, size_t col_idx) {

    size_t length = vec.size();

    assert(mat[0].size() > col_idx);
    assert(mat.size() == length);

    for (size_t i = 0; i < length; i++) {
        mat[i][col_idx] = vec.get(i);
    }
}

Vector
getCol(const vector<vector<double>>& mat, size_t col_idx) {
    size_t dim = mat.size();

    assert(mat[0].size() > col_idx);

    Vector col(mat.size());

    for (size_t i = 0; i < dim; i++) {
        col.set(i, mat[i][col_idx]);
    }

    return col;
}

void
printMatrix(vector<vector<double>>& mat,
            size_t row_start, size_t row_end,
            size_t col_start, size_t col_end) {
    size_t m = mat.size();
    size_t n = mat[0].size();

    assert(row_start >= 0);
    assert(col_start >= 0);
    assert(row_end <= m);
    assert(col_end <= n);

    char buf[1024];
    sprintf(buf, "");

    for (int i = row_start; i < row_end; i++) {
        for (int j = col_start; j < col_end; j++) {
            if (j == col_start)
                sprintf(buf, "%s[", buf);

            sprintf(buf, "%s%.4f", buf, mat[i][j]);

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
leastSquare(vector<vector<double>>& H, size_t size, double beta) {
    Eigen::MatrixXf A(size+1, size);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(size+1);
    Eigen::VectorXf y_(size);
    Vector y(size);

    b(0) = beta;

    for (int i = 0; i < size+1; i++) {
        for (int j = 0; j < size; j++) {
            A(i, j) = H[i][j];
        }
    }


    y_ = (A.transpose() * A).ldlt().solve(A.transpose() * b);

    for (int i = 0; i < size; i++) {
        y.set(i, y_(i));
    }

    return y;
}

Vector
gmres(const vector<vector<double>>& A,
      const Vector& b,
      size_t m, size_t tol, size_t maxit) {

    size_t dim = A[0].size();
    size_t nit = 0;
    Vector x0(dim);

    assert(dim == b.size());

    // Use trivial preconditioner for now
    auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    while (nit < maxit) {
        auto H = generateMatrix(m+1, m);
        auto Z = generateMatrix(dim, m);
        auto V = generateMatrix(dim, m+1);

        // TODO: the GMRES algorithm
        auto r0 = b.sub(mvAx(A, x0));
        double beta = r0.norm2();
        Vector x(dim);
        setCol(V, r0.mulS(1.0 / beta), 0);

        // TODO: get krylov space
        for(int j = 0; j < m; j++) {
            setCol(Z, mvAx(Prcnd, getCol(V, j)), j);
            auto w = mvAx(A, getCol(Z, j));

            for (int i = 0; i < j; i++) {
                auto v = getCol(V, i);
                H[i][j] = w.dotV(v);
                w.isub(v.mulS(H[i][j]));
            }

            H[j+1][j] = w.norm2();
            setCol(V, w.mulS(1.0 / H[j+1][j]), j+1);

            // printMatrix(H, 0, j+1, 0, j);
            // cout << "---" << endl;

            auto y = leastSquare(H, j+1, beta);
            auto x = x0.add(mvZy(Z, y, j+1));

            auto res_norm = mvAx(A, x).sub(b).norm2();

            cout << "#" << j <<
                " Residual=" << res_norm << endl;

        }

        x0 = x;
        nit++;

        break;
    }
    return x0;
}

int main(int argc, char *argv[])
{
    cout << "A = ../data/cage4.mtx" << endl;

    auto A = loadMTXFile("../data/cage4.mtx");

    auto b = Vector(A.size());
    b.set(0, 1.0);

    gmres(A, b, 10, 0, 10);

    return 0;
}
