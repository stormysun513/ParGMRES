#include <iostream>
#include <algorithm>
#include <cassert>

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <cusp/csr_matrix.h>
#include <cusp/io/matrix_market.h>
#include <cusp/krylov/gmres.h>
#include <cusp/print.h>

#include "gmres.cuh"
#include "CycleTimer.h"

#include "utils.h"

#define WARP_SIZE               32
#define MAX_THREAD_NUM          1024
#define KRYLOV_M                100
#define INDEX( i, j, dim )      (i * dim + j)
#define ROUND( num, base )      ((num + base - 1)/base)
#define HANDLE_ERROR( err )     (cuda_handle_error(err, __FILE__, __LINE__))

/*
__inline__ __device__
float warp_reduce_sum(float val) {

    for (int offset = WARP_SIZE/2; offset > 0; offset /= 2)
        val += __shfl_down(val, offset);
    return val;
}

__global__
void device_reduce_warp_atomic_kernel(float *in, float* out, int N) {

    float sum = float(0);
    for(int i = blockIdx.x * blockDim.x + threadIdx.x;
            i < N;
            i += blockDim.x * gridDim.x) {
        sum += in[i];
    }
    sum = warp_reduce_sum(sum);
    if (threadIdx.x & (WARP_SIZE - 1) == 0)
        atomicAdd(out, sum);
}
*/

/* kernel functions */

__global__
void s_mem_init(float *addr, float value, int N) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i >= N) return;
    addr[i] = value;
}

__global__
void s_x_sqrt(float *res, float *x, int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;
    res[i] = sqrt(x[i]);
}

__global__
void s_x_div_a(float *res, float *x, float *a, int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;
    res[i] = x[i] / (*a);
}

__global__
void s_x_dot_y(float *res, float *x, float *y, int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;
    float value = x[i] * y[i];
    atomicAdd(res, value);
}

__global__
void s_x_sub_ay(float *x, float *y, float *a, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;
    x[i] = x[i] - y[i] * (*a);
}

__global__
void gmres_update_x(float *x, float *V, float *y, int m, int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float entry = .0;

    if(i >= N) return;

    for(int k = 0; k < m; k++){
        entry += V[k*N+i] * y[k];
    }
    x[i] += entry;
}

__global__
void s_mat_mul_x(float *res, csr_mat_t mat, float *x){

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    int nrow = mat.nrow;
    int *rowstart = mat.rowstart;
    int *cindex = mat.cindex;
    float *value = mat.value;

    if(i >= nrow) return;

    int start_idx = rowstart[i];
    int end_idx = rowstart[i+1];

    float temp = 0.0;
    for (int k = start_idx; k < end_idx; ++k) {
        int j = cindex[k];
        temp += value[k] * x[j];
    }
    res[i] = temp;
}

__global__
void gmres_compute_r0(float *r0, csr_mat_t mat, float *x, vec_t vec, float *beta){

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    int nrow = mat.nrow;
    int *rowstart = mat.rowstart;
    int *cindex = mat.cindex;
    float *value = mat.value;

    if(i >= nrow) return;

    int start_idx = rowstart[i];
    int end_idx = rowstart[i+1];

    float temp = 0.0;
    for (int k = start_idx; k < end_idx; ++k) {
        int j = cindex[k];
        temp += value[k] * x[j];
    }
    r0[i] = temp - vec.value[i];

    float square = r0[i];
    square *= square;
    atomicAdd(beta, square);
}

/* host functions */

void cuda_handle_error(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}

void mem_log(float* device, int N) {

    float host[1024];
    char buf[4096];
    int min = std::min(1024, N);

    assert(min >= 0);

    HANDLE_ERROR(cudaMemcpy(host, device, N*sizeof(float), cudaMemcpyDeviceToHost));

    for(int i = 0; i < min; i++){
        sprintf(buf, "%s%lf, ", buf, host[i]);
    }
    std::cout << buf << std::endl;
}

static void display_gpu_info() {

    const int kb = 1024;
    const int mb = kb * kb;

    std::cout << "\nCUDA version: v" << CUDART_VERSION << std::endl;
    std::cout << "Thrust version: v" << THRUST_MAJOR_VERSION << ".";
    std::cout << THRUST_MINOR_VERSION << std::endl;

    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));

    std::cout << "\nCUDA Devices: \n\n";

    for(int i = 0; i < devCount; ++i) {
        cudaDeviceProp props;
        HANDLE_ERROR(cudaGetDeviceProperties(&props, i));
        std::cout << i << ":\n  " << props.name << ": " << props.major << "." << props.minor << std::endl;
        std::cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << std::endl;
        std::cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << std::endl;
        std::cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << std::endl;
        std::cout << "  Block registers: " << props.regsPerBlock << std::endl << std::endl;

        std::cout << "  Warp size:         " << props.warpSize << std::endl;
        std::cout << "  Threads per block: " << props.maxThreadsPerBlock << std::endl;
        std::cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", "
                                             << props.maxThreadsDim[1] << ", "
                                             << props.maxThreadsDim[2] << " ]" << std::endl;
        std::cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", "
                                             << props.maxGridSize[1] << ", "
                                             << props.maxGridSize[2] << " ]" << std::endl;
        std::cout << std::endl;
    }
}

void gmres(csr_mat_t mat, vec_t vec, int m, float tol, int maxit){

    int dim = mat.ncol;
    int nit = 0;
    int innit = 0;
    int outnit = 0;

    float *H;
    float *V;

    float *x;
    float *y;
    float *w;
    float *r0;

    float *beta;
    float *tmp1;
    float *tmp2;

    size_t H_bytes = (m+1) * m * sizeof(float);
    size_t y_bytes = m * sizeof(float);
    float *H_host_data = (float *)calloc((m+1) * m, sizeof(float));
    float *y_host_data = (float *)calloc(m, sizeof(float));
    float *beta_host = (float *)malloc(sizeof(float));

    std::cout << "\nOur GMRES solution:\n\n";

    HANDLE_ERROR(cudaMalloc((void**)&H, (m+1) * m * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&V, (m+1) * dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&x, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&y, m * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&w, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&r0, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&beta, sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&tmp1, sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&tmp2, sizeof(float)));

    int blocks = ROUND(dim, MAX_THREAD_NUM);
    int threads = MAX_THREAD_NUM;

    s_mem_init<<<blocks, threads>>>(x, .0, dim);

    while(nit < maxit){

        // kernel 1: compute r0 and beta
        gmres_compute_r0<<<blocks, threads>>>(r0, mat, x, vec, tmp1);
        s_x_sqrt<<<1,1>>>(beta, tmp1, 1);
        s_x_div_a<<<blocks, threads>>>(V, r0, beta, dim);

        innit = 0;

        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {

            // tStart = CycleTimer::currentSeconds();

            // compute mat and vec mulplication (mat, V, j),  w can be placed at V(:, j+1) in the future
            s_mat_mul_x<<<blocks, threads>>>(w, mat, (V + j*dim));

            // for (size_t i = 0; i < j; i++) {
            //     Vector v = V.getCol(i);
            //     H.set(i, j, w.dotV(v));
            //     w.isub(v.mulS(H.get(i, j)));
            // }
            for(size_t i = 0; i < j; i++){
                s_x_dot_y<<<blocks, threads>>>(H+i*(KRYLOV_M+1)+j, w, (V+i*dim), dim);
                s_x_sub_ay<<<blocks, threads>>>(w, (V+i*dim), H+i*(KRYLOV_M+1)+j, dim);
            }

            //H.set(j+1, j, w.norm2());
            //V.setCol(j+1, w.mulS(1.0 / H.get(j+1, j)));
            float *out = H+(j+1)*(KRYLOV_M+1)+j;
            s_x_dot_y<<<blocks, threads>>>(out, w, w, dim);
            s_x_sqrt<<<1,1>>>(tmp1, out, 1);
            s_x_div_a<<<blocks, threads>>>(V+(j+1)*dim, w, tmp1, dim);


            // tKrylov += CycleTimer::currentSeconds() - tStart;
            // tStart = CycleTimer::currentSeconds();

	        // TODO: make cuda version LLS
            // Vector y = leastSquareWithQR(H, j+1, beta);
            cudaMemcpy(H_host_data, H, H_bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(beta_host, beta, sizeof(float), cudaMemcpyDeviceToHost);
            Matrix H_host_mat(m+1, m, H_host_data);
            Vector y_host_vec = leastSquareWithQR(H_host_mat, j+1, *beta_host);
            cudaMemcpy(y, y_host_vec.getData(), y_bytes, cudaMemcpyHostToDevice);

            // x = x0.add(V.mulPartialT(y, j+1));
            // float res_norm = A.mul(x).sub(b).norm2();
            gmres_update_x<<<blocks, threads>>>(x, V, y, j+1, dim);
            gmres_compute_r0<<<blocks, threads>>>(r0, mat, x, vec, tmp1);
            s_x_sqrt<<<1,1>>>(tmp2, tmp1, 1);

            nit++;
            innit++;

            // tLLS += CycleTimer::currentSeconds() - tStart;

            // if (res_norm < tol * b.norm2()) {
            //     cout << "FGMRES converged to relative tolerance: "
            //          << res_norm / b.norm2() << " at iteration " << nit
            //          << " (out: " << outnit << ", in: " << innit << ")" << endl;

            //     sprintf(buf, "[%.3f] ms in Krylov \n", tKrylov * 1000);
            //     cout << buf;
            //     sprintf(buf, "[%.3f] ms in LLS \n", tLLS * 1000);
            //     cout << buf;

            //     return x;
            // }
        }
        outnit++;
    }

    HANDLE_ERROR(cudaFree(H));
    HANDLE_ERROR(cudaFree(V));
    HANDLE_ERROR(cudaFree(x));
    HANDLE_ERROR(cudaFree(y));
    HANDLE_ERROR(cudaFree(w));
    HANDLE_ERROR(cudaFree(r0));
    HANDLE_ERROR(cudaFree(tmp1));
    HANDLE_ERROR(cudaFree(beta));
}

template <typename Matrix, typename Vector>
static void gmres_ref(const Matrix& A, Vector& x, const Vector& b){

    float start_time;
    float end_time;
    char buf[1024];

    // reference answer
    std::cout << "\nReference answer from CUSP library:\n\n";

    // set stopping criteria:
    cusp::monitor<float> monitor(b, KRYLOV_M, 1e-6, 0, false);
    int restart = 50;

    // run gmres to solve
    start_time = CycleTimer::currentSeconds();
    cusp::krylov::gmres(A, x, b, restart, monitor);
    end_time = CycleTimer::currentSeconds();

    // print the performance
    sprintf(buf, "[%.3f] ms in total (CSR Sparse GMRES) \n\n",
            (end_time - start_time) * 1000);
    std::cout << buf;

    // print out result for debug
    cusp::print(x);
}

void run(void) {

    display_gpu_info();

    // create an empty sparse matrix structure (CSR format)
    cusp::csr_matrix<int, float, cusp::device_memory> A;

    // load a matrix stored in MatrixMarket format
    cusp::io::read_matrix_market_file(A, "../../data/cage4.mtx");

    // allocate storage for solution (x) and right hand side (b)
    cusp::array1d<float, cusp::device_memory> x(A.num_cols);     // 0
    cusp::array1d<float, cusp::device_memory> b(A.num_rows, 1);     // 1

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

    gmres_ref(A, x, b);
    gmres(csr_mat, vec, 100, 1e-6, 1000);
}
