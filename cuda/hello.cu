#include <iostream>                                                                                                                                                                 
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <cuda.h>

using namespace std;

int main(){
    size_t N = 10; 
    int * raw_ptr;

    cudaMalloc((void **) &raw_ptr, N * sizeof(int));
    thrust::device_ptr<int> dev_ptr = thrust::device_pointer_cast(raw_ptr);
    thrust::fill(dev_ptr, dev_ptr + N, (int) 0); 

    cudaFree(raw_ptr);

    cout << "finished" << endl; 

    return 0;
}
