#include "main.hpp"

#include <iostream>
#include <vector>
#include <cassert>

#define MAX_KRYLOV_DIM  500
#define MAX_ITERS       100

using namespace std;

template <class T>
vector<T>
mvAx(const vector<vector<T>>& A, const vector<T>& x){
    vector<T> b(x.size());

    for (size_t i = 0; i < A.size(); ++i) {
        b[i] = 0;
        for (size_t j = 0; j < A[0].size(); ++j) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}

template <class T>
vector<T>
vSub(const vector<T>& a, const vector<T>& b){
    vector<T> ans(a.size());

    for(size_t i = 0; i < a.size(); i++){
        ans[i] = a[i] - b[i];
    }
    return ans;
}

vector<double>
gmres(const vector<vector<double>>& A, 
        const vector<double>& b, 
        size_t m, 
        size_t tol, 
        size_t maxit){

    vector<double> ans(A[0].size());
    size_t nit = 0;

    m = (m > MAX_KRYLOV_DIM) ? MAX_KRYLOV_DIM : m;
    maxit = (maxit > MAX_ITERS) ? MAX_ITERS : maxit;

    assert(m > 0);
    assert(maxit > 0);

    while(nit < maxit){
        vector<vector<double>> H(m+1);
        for(size_t i = 0; i < m+1; i++){
            H[i] = vector<double>(m);
        }
        // TODO: the GMRES algorithm       

        nit++;
    }
    return ans;
} 

int main(int argc, char *argv[])
{
    cout << "Hello, world!" << endl;
    return 0;
}
