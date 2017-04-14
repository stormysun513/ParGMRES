#include "main.hpp"
#include "utils.hpp"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100

using namespace std;

Vector
mvAx(const vector<vector<double>>& A, const Vector& x) {
    Vector b(x.size());

    for (size_t i = 0; i < A.size(); ++i) {
        b.set(i, 0);
        for (size_t j = 0; j < A[0].size(); ++j) {
            b.set(i, b.get(i) + A[i][j] * x.get(j));
        }
    }
    return b;
}

Vector
vMul(const Vector& vec, double scale) {
    Vector ans(vec.size());

    for (size_t i = 0; i < vec.size(); i++){
        ans.set(i, scale * vec.get(i));
    }
    return ans;
}

Vector
vSub(const Vector& vec1, const Vector& vec2) {
    Vector ans(vec1.size());

    for (size_t i = 0; i < vec1.size(); i++) {
        ans.set(i, vec1.get(i) - vec2.get(i));
    }
    return ans;
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
getCol(vector<vector<double>>& mat, size_t col_idx) {
    size_t dim = mat.size();

    assert(dim > col_idx);

    Vector col(mat.size());

    for (size_t i = 0; i < dim; i++) {
        col.set(i, mat[i][col_idx]);
    }

    return col;
}

Vector
gmres(const vector<vector<double>>& A,
      const Vector& b,
      size_t m, size_t tol, size_t maxit) {

    size_t dim = A[0].size();
    size_t nit = 0;
    Vector x0(dim);

    // Use trivial preconditioner for now
    auto Prcnd = identityMatrix(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    while (nit < maxit) {
        auto H = generateMatrix(m+1, m);
        auto Z = generateMatrix(m, dim);
        auto V = generateMatrix(m+1, dim);

        // TODO: the GMRES algorithm
        auto r0 = vSub(b, mvAx(A, x0));
        double beta = r0.norm2();
        Vector x(dim);
        setCol(V, vMul(r0, 1.0 / beta), 0);

        // TODO: get krylov space
        for(int j = 0; j < m; j++){
            setCol(Z, mvAx(Prcnd, getCol(V, j)), j);
        }

        x0 = x;
        nit++;
    }
    return x0;
}

int main(int argc, char *argv[])
{
    cout << "Hello, world!" << endl;
    return 0;
}
