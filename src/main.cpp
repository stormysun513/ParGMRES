#include "main.hpp"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100

using namespace std;

vector<double>
mvAx(const vector<vector<double>>& A, const vector<double>& x){
    vector<double> b(x.size());

    for (size_t i = 0; i < A.size(); ++i) {
        b[i] = 0;
        for (size_t j = 0; j < A[0].size(); ++j) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}

vector<double>
vMul(const vector<double>& vec, double scale){
    vector<double> ans(vec.size());

    for(size_t i = 0; i < vec.size(); i++){
        ans[i] = scale * vec[i];
    }
    return ans;
}

vector<double>
vSub(const vector<double>& vec1, const vector<double>& vec2){
    vector<double> ans(vec1.size());

    for(size_t i = 0; i < vec1.size(); i++){
        ans[i] = vec1[i] - vec2[i];
    }
    return ans;
}

double norm(const vector<double>& vec){
    double res = .0f;

    for(double num: vec){
        res += num;
    }
    return sqrt(res);
}

vector<vector<double>>
generateMatrix(int m, int n){
    vector<vector<double>> matrix(m);

    for(size_t i = 0; i < m; i++){
        matrix[i] = vector<double>(n);
    }
    return matrix;
}

void setRow(vector<vector<double>>& mat, const vector<double>& vec, int row){

    size_t length = vec.size();

    assert(mat.size() > row);
    assert(mat[0].size() == length);

    for(int i = 0; i < length; i++){
        mat[row][i] = vec[i];
    }
}

vector<double>
gmres(const vector<vector<double>>& A, 
        const vector<double>& b, 
        size_t m, 
        size_t tol, 
        size_t maxit){

    size_t dim = A[0].size();
    size_t nit = 0;
    vector<double> x0(dim);

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    while(nit < maxit){
        auto H = generateMatrix(m+1, m);
        auto Z = generateMatrix(dim, m);
        auto V = generateMatrix(m+1, dim);
         
        // TODO: the GMRES algorithm       
        auto r0 = vSub(b, mvAx(A, x0));
        double beta = norm(r0);
        vector<double> x(dim);
        setRow(V, vMul(r0, 1/beta), 0);
 
        // TODO: get krylov space
        for(int j = 0; j < m; j++){
            
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
