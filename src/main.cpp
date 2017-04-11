#include "main.hpp"

#include <iostream>
#include <vector>

using namespace std;

vector<double>
mvAx(vector<vector<double>> A, vector<double> x) {
    vector<double> b(x.size());

    for (auto i = 0; i < A.size(); ++i) {
        for (auto j = 0; j < A[0].size(); ++j) {
            b[i] += A[i][j] * x[j];
        }
    }

    return b;
}

int main(int argc, char *argv[])
{
    cout << "Hello, world!" << endl;
    return 0;
}
