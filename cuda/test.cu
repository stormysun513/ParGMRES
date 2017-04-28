#include <iostream>

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <cusp/csr_matrix.h>
#include <cusp/io/matrix_market.h>
#include <cusp/krylov/gmres.h>
#include <cusp/print.h>

#include "CycleTimer.h"

using namespace std;

#define MAX_THREAD_PER_BLOCK    1024 

typedef struct csr_mat_t {
    int *rowstart;
    int *cindex;
    double *value;
    int nrow;
    int ncol;
    int nnz;
} csr_mat_t;

typedef struct vec_t {
    double *value;
    int size;
} vec_t;

__global__ void init_b_vec(vec_t vec){
    
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    int size = vec.size;
    double *b = vec.value;

    if(idx < size){
        b[idx] = idx;
    }
}

void display_gpu_info() {
    
    const int kb = 1024;
    const int mb = kb * kb;
    
    cout << endl << endl;
    cout << "CUDA version:   v" << CUDART_VERSION << endl;
    cout << "Thrust version: v" << THRUST_MAJOR_VERSION << ".";
    cout << THRUST_MINOR_VERSION << endl << endl; 

    int devCount;
    cudaGetDeviceCount(&devCount);
    
    cout << "CUDA Devices: " << endl << endl;

    for(int i = 0; i < devCount; ++i) {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        cout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
        cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
        cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
        cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
        cout << "  Block registers: " << props.regsPerBlock << endl << endl;

        cout << "  Warp size:         " << props.warpSize << endl;
        cout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
        cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " 
                                             << props.maxThreadsDim[1] << ", " 
                                             << props.maxThreadsDim[2] << " ]" << endl;
        cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " 
                                             << props.maxGridSize[1] << ", " 
                                             << props.maxGridSize[2] << " ]" << endl;
        cout << endl;
    }
}

void gmres(csr_mat_t mat, vec_t vec, int m, double tol, int maxit){
    
}

int main(void)
{
    double start_time;
    double end_time;
    char buf[1024];

    display_gpu_info();

    // create an empty sparse matrix structure (CSR format)
    cusp::csr_matrix<int, double, cusp::device_memory> A;

    // load a matrix stored in MatrixMarket format
    cusp::io::read_matrix_market_file(A, "../data/cage4.mtx");

    // allocate storage for solution (x) and right hand side (b)
    cusp::array1d<double, cusp::device_memory> x(A.num_cols, 0);
    cusp::array1d<double, cusp::device_memory> b(A.num_rows, 1);

    // get raw pointer
    csr_mat_t csr_mat;
    csr_mat.nrow = A.num_rows;
    csr_mat.ncol = A.num_cols;
    csr_mat.nnz = A.values.size();
    csr_mat.value = thrust::raw_pointer_cast(A.values.data());
    csr_mat.cindex = thrust::raw_pointer_cast(A.column_indices.data());
    csr_mat.rowstart = thrust::raw_pointer_cast(A.row_offsets.data());

    vec_t vec;
    vec.value = thrust::raw_pointer_cast(b.data());
    vec.size = b.size();
    
    // initialize b vector
    dim3 gridDims(1);
    init_b_vec<<<gridDims, vec.size>>>(vec);

    // set stopping criteria:
    cusp::monitor<double> monitor(b, 100, 1e-6, 0, false);
    int restart = 50;

    // run gmres to solve 
    start_time = CycleTimer::currentSeconds();
    cusp::krylov::gmres(A, x, b, restart, monitor);
    end_time = CycleTimer::currentSeconds();

    // print the performance
    sprintf(buf, "[%.3f] ms in total (CSR Sparse GMRES) \n\n", 
            (end_time - start_time) * 1000);
    cout << buf;

    // print out result for debug
    cusp::print(A);
    cusp::print(b);
    cusp::print(x);

    return 0;
}


