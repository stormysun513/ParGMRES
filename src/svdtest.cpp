#include <iostream>

#include "CycleTimer.h"

#include "utils.h"

int main(){
    
    // Test 1
    // Matrix A = loadMTXToMatrix("../data/cage4.mtx");
    // Matrix U, V;
    // Vector sigma;

    // svdJacobiMethod(U, V, sigma, A);

    // std::cout << "\nSigma: " << std::endl;
    // printVector(sigma);
    // std::cout << std::endl;

    // std::cout << "\n\nV matrix: " << std::endl;
   
    // Test 2 
    //Matrix H = loadMTXToMatrix("../data/H.mtx");
    //size_t j = 10;
    //double beta = 1.82896;

    //Vector y1 = leastSquareWithEigen(H, j+1, beta);
    //Vector y2 = leastSquareWithPowerMethod(H, j+1, beta);
    //Vector y3 = leastSquareWithJacobi(H, j+1, beta);

    //printVector(y1);
    //printVector(y2);
    //printVector(y3);

    // Test 3
    char buf[1024];
    double start_time;
    double end_time;

    //for(size_t i = 10000; i <= 10000000; i *= 10){

    //    Vector a = randUnitUniformVector(i); 

    //    std::cout << "L2 Norm test: " << i << std::endl;

    //    start_time = CycleTimer::currentSeconds();
    //    a.norm2();
    //    end_time = CycleTimer::currentSeconds();

    //    sprintf(buf, "[%.3f] ms (Sequential) \n", 
    //            (end_time - start_time) * 1000);
    //    std::cout << buf;

    //    start_time = CycleTimer::currentSeconds();
    //    l2norm(a);
    //    end_time = CycleTimer::currentSeconds();

    //    sprintf(buf, "[%.3f] ms (Omp) \n\n", 
    //            (end_time - start_time) * 1000);
    //    std::cout << buf;
    //}


    for(size_t i = 100; i <= 1000; i += 100){

        Vector y = randUnitUniformVector(i); 
        Vector x = randUnitUniformVector(i); 
        Matrix A = randUniformMatrix(i); 

        std::cout << "mvmul test: " << i << std::endl;

        start_time = CycleTimer::currentSeconds();
        A.mul(x);
        end_time = CycleTimer::currentSeconds();

        sprintf(buf, "[%.3f] ms (Sequential) \n", 
                (end_time - start_time) * 1000);
        std::cout << buf;

        start_time = CycleTimer::currentSeconds();
        matVecMul(y, A, x);
        end_time = CycleTimer::currentSeconds();

        sprintf(buf, "[%.3f] ms (Omp) \n\n", 
                (end_time - start_time) * 1000);
        std::cout << buf;
    }

    return 0;
}
