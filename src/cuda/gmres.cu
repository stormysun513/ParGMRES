#include <iostream>
#include <algorithm>
#include <cmath>
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
#define MAX_THREAD_NUM          256
#define KRYLOV_M                100
#define INDEX( i, j, dim )      (i * dim + j)
#define ROUND( num, base )      ((num + base - 1)/base)
#define HANDLE_ERROR( err )     (cuda_handle_error(err, __FILE__, __LINE__))

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
    printf("i = %d, x[i] = %f \n", i, x[i]);
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
    printf("value = %f\n", value);
    printf("res = %f\n", res);
}

__global__
void s_x_sub_ay(float *x, float *y, float *a, int N){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= N) return;
    x[i] = x[i] - y[i] * (*a);
}

__global__
void gmres_update_x(float *x, float *x0, float *V, float *y, int m, int N){

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    float entry = .0;

    if(i >= N) return;

    for(int k = 0; k < m; k++){
        entry += V[k*N+i] * y[k];
    }
    x[i] = x0[i] + entry;
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
void gmres_compute_r0(float *r0, csr_mat_t mat, float *x, vec_t vec){

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
    r0[i] = vec.value[i] - temp;
}

/* thrust related subroutine */

template <class T>
struct square
{
    __host__ __device__
    T operator()(const T& x) const {
        return x * x;
    }
};

template <class T>
struct saxpy_functor: public thrust::binary_function<T, T, T>
{
    const T a;
    saxpy_functor(T _a) : a(_a) {}

    __host__ __device__
    T operator()(const T& x, const T& y) const {
        return a * x + y;
    }
};

template <class T>
struct sysax_functor: public thrust::binary_function<T, T, T>
{
    const T a;
    sysax_functor(T _a) : a(_a) {}

    __host__ __device__
    T operator()(const T& x, const T& y) const {
        return y - a * x;
    }
};

template <class T>
struct sax_functor: public thrust::unary_function<T, T> // Arguement, Result
{
    const T a;
    sax_functor(T _a) : a(_a) {}

    __host__ __device__
    T operator()(const T& x) const {
        return a * x;
    }
};

template <class Iteratable>
void saxpy(float a, Iteratable& X, Iteratable& Y) {
    // Y <- A * X + Y
    thrust::transform(X.begin(), X.end(), Y.begin(), Y.begin(), saxpy_functor<float>(a));
}

template <class Iterator>
void sax(float a, Iterator xbegin, Iterator xend, Iterator ybegin) {
    // Y <- A * X
    thrust::transform(xbegin, xend, ybegin, sax_functor<float>(a));
}

template <class Iterator>
void sysax(float a, Iterator xbegin, Iterator xend, Iterator ybegin) {
    // Y <- Y - A * X
    thrust::transform(xbegin, xend, ybegin, ybegin, sysax_functor<float>(a));
}

template <class Iterator>
static float l2norm(Iterator begin, Iterator end){

    square<float>        unary_op;
    thrust::plus<float> binary_op;
    float init = 0;

    return std::sqrt(thrust::transform_reduce(begin, end, unary_op, init, binary_op));
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

    buf[0] = '\0';

    HANDLE_ERROR(cudaMemcpy(host, device, N*sizeof(float), cudaMemcpyDeviceToHost));

    for(int i = 0; i < min; i++){
        sprintf(buf, "%s%lf, ", buf, host[i]);
    }
    std::cout << buf << std::endl;
}

void mem_log_2d(float* device, int M, int N, int dim) {

    float host[2048];
    char buf[4096];
    int min = std::min(1024, N);

    assert(min >= 0);

    buf[0] = '\0';

    HANDLE_ERROR(cudaMemcpy(host, device, M*dim*sizeof(float), cudaMemcpyDeviceToHost));

    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            sprintf(buf, "%s%lf, ", buf, host[i*dim + j]);
        }
        sprintf(buf, "%s\n", buf);
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

void gmres(csr_mat_t mat, vec_t res, vec_t vec, int m, float tol, int maxit){

    int dim = mat.ncol;
    int nit = 0;
    int innit = 0;
    int outnit = 0;

    bool terminate = false;

    float *H;
    float *V;

    float *x;
    float *y;
    float *w;
    float *r0;
    float *x0;

    float *beta;
    float *tmp1;
    float *tmp2;

    size_t H_bytes = (m+1) * m * sizeof(float);
    size_t y_bytes = m * sizeof(float);
    float *H_host_data = (float *)calloc((m+1) * m, sizeof(float));
    float *y_host_data = (float *)calloc(m, sizeof(float));
    float *beta_host = (float *)malloc(sizeof(float));

    float *temp_host_float = (float *)malloc(sizeof(float));

    std::cout << "\nOur GMRES solution:\n\n";

    HANDLE_ERROR(cudaMalloc((void**)&H, (m+1) * m * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&V, (m+1) * dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&x, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&x0, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&y, m * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&w, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&r0, dim * sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&beta, sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&tmp1, sizeof(float)));
    HANDLE_ERROR(cudaMalloc((void**)&tmp2, sizeof(float)));

    int blocks = ROUND(dim, MAX_THREAD_NUM);
    int threads = MAX_THREAD_NUM;

    thrust::device_ptr<float> dp_b = thrust::device_pointer_cast(vec.value);
    thrust::device_ptr<float> dp_x = thrust::device_pointer_cast(x);
    thrust::device_ptr<float> dp_x0 = thrust::device_pointer_cast(x0);
    thrust::device_ptr<float> dp_tmp1 = thrust::device_pointer_cast(tmp1);
    thrust::device_ptr<float> dp_H = thrust::device_pointer_cast(H);

    float b_norm2 = l2norm(dp_b, dp_b + dim);

    thrust::fill(dp_x0, dp_x0+dim, .0);
    thrust::fill(dp_H, dp_H + (m+1)*m, .0);

    while(nit < maxit) {

        // kernel 1: compute r0 and beta
        gmres_compute_r0<<<blocks, threads>>>(r0, mat, x0, vec);

        // s_x_sqrt<<<1,1>>>(beta, tmp1, 1);
        thrust::device_ptr<float> dp_r0 = thrust::device_pointer_cast(r0);
        *temp_host_float = l2norm(dp_r0, dp_r0 + dim);

        cudaMemcpy(beta, temp_host_float, sizeof(float), cudaMemcpyHostToDevice);

        s_x_div_a<<<blocks, threads>>>(V, r0, beta, dim);

        //float res = std::sqrt(thrust::reduce(dp_tmp1, dp_tmp1+dim));

        innit = 0;

        // Generate krylov subspace
        for(size_t j = 0; j < m; j++) {

            // tStart = CycleTimer::currentSeconds();

            // compute mat and vec mulplication (mat, V, j),  w can be placed at V(:, j+1) in the future
            s_mat_mul_x<<<blocks, threads>>>(w, mat, (V + j*dim));
            cudaDeviceSynchronize();

            // DEBUG
            if(nit < 3) {
                std::cout << "w:\n";
                mem_log(w, dim);

                std::cout << "V:\n";
                mem_log_2d(V, j+1, dim, dim);
            }

            // for (size_t i = 0; i < j; i++) {
            //     Vector v = V.getCol(i);
            //     H.set(i, j, w.dotV(v));
            //     w.isub(v.mulS(H.get(i, j)));
            // }
            for(size_t i = 0; i < j; i++){

                // s_x_dot_y<<<blocks, threads>>>(H+i*(KRYLOV_M+1)+j, w, (V+i*dim), dim);
                thrust::device_ptr<float> dp_w = thrust::device_pointer_cast(w);
                thrust::device_ptr<float> dp_vi = thrust::device_pointer_cast(V+i*dim);
                float *dot = (float *)malloc(sizeof(float));
                *dot = thrust::inner_product(dp_w, dp_w+dim, dp_vi, 0.f);
                float *hij = H+i*m+j;

                cudaMemcpy(hij, dot, sizeof(float), cudaMemcpyHostToDevice);
                s_x_sub_ay<<<blocks, threads>>>(w, (V+i*dim), hij, dim);
                cudaDeviceSynchronize();

                free(dot);
            }

            //H.set(j+1, j, w.norm2());
            //V.setCol(j+1, w.mulS(1.0 / H.get(j+1, j)));
            thrust::device_ptr<float> dp_w = thrust::device_pointer_cast(w);

            // s_x_dot_y<<<blocks, threads>>>(out, w, w, dim);
            // s_x_sqrt<<<1,1>>>(tmp1, out, 1);
            *temp_host_float = l2norm(dp_w, dp_w + dim);
            cudaMemcpy(tmp1, temp_host_float, sizeof(float), cudaMemcpyHostToDevice);

            cudaMemcpy(H+(j+1)*m+j, temp_host_float, sizeof(float), cudaMemcpyHostToDevice);
            s_x_div_a<<<blocks, threads>>>(V+(j+1)*dim, w, tmp1, dim);


            // tKrylov += CycleTimer::currentSeconds() - tStart;
            // tStart = CycleTimer::currentSeconds();

            if(nit < 3){
                std::cout << "H:\n";
                mem_log_2d(H, j+2, j+1, m);

                std::cout << "V:\n";
                mem_log_2d(V, j+2, dim, dim);
            }


	        // TODO: make cuda version LLS
            // Vector y = leastSquareWithQR(H, j+1, beta);
            cudaMemcpy(H_host_data, H, H_bytes, cudaMemcpyDeviceToHost);
            cudaMemcpy(beta_host, beta, sizeof(float), cudaMemcpyDeviceToHost);
            Matrix H_host_mat(m+1, m, H_host_data);
            Vector y_host_vec = leastSquareWithQR(H_host_mat, j+1, *beta_host);
            cudaMemcpy(y, y_host_vec.getData(), y_bytes, cudaMemcpyHostToDevice);

            if(nit < 3){
                std::cout << "y:\n";
                printVector(y_host_vec);
            }

            // x = x0.add(V.mulPartialT(y, j+1));
            // float res_norm = A.mul(x).sub(b).norm2();
            gmres_update_x<<<blocks, threads>>>(x, x0, V, y, j+1, dim);


            if(nit < 3){
                std::cout << "x:\n";
                mem_log(x, dim);
            }


            gmres_compute_r0<<<blocks, threads>>>(r0, mat, x, vec);
            float res_norm = l2norm(dp_r0, dp_r0 + dim);

            nit++;
            innit++;

            // tLLS += CycleTimer::currentSeconds() - tStart;

            if (res_norm < tol * b_norm2) {
                std::cout << "FGMRES converged to relative tolerance: "
                     << res_norm / b_norm2 << " at iteration " << nit
                     << " (out: " << outnit << ", in: " << innit << ")" << std::endl;

                cudaMemcpy(res.value, x, dim * sizeof(float), cudaMemcpyDeviceToHost);

                // sprintf(buf, "[%.3f] ms in Krylov \n", tKrylov * 1000);
                // std::cout << buf;
                // sprintf(buf, "[%.3f] ms in LLS \n", tLLS * 1000);
                // std::cout << buf;

                terminate = true;
                break;
            }
        }

        if (terminate)
            break;

        // x0 = x;
        cudaMemcpy(x0, x, dim * sizeof(float), cudaMemcpyDeviceToDevice);
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

template <typename Mat, typename Vec>
static void gmres_new(Mat& A, Vec& x, Vec& b){

    float t_start;
    float t_temp;
    float t_krylov = 0;
    float t_lls = 0;
    float t_overall = 0;
    char buf[1024];

    int m = 100;
    int maxit = 10000;
    float tol = 1e-6;

    int nit = 0;
    int innit = 0;
    int outnit = 0;
    bool terminate = false;

    int nrow = A.num_rows;
    int ncol = A.num_cols;

    int blocks = ROUND(nrow, MAX_THREAD_NUM);
    int threads = MAX_THREAD_NUM;

    size_t H_bytes = (m+1) * m * sizeof(float);
    size_t y_bytes = m * sizeof(float);
    float *H_host_data = (float *)calloc((m+1) * m, sizeof(float));
    float *y_host_data = (float *)calloc(m, sizeof(float));
    float *beta_host = (float *)malloc(sizeof(float));

    thrust::device_vector<float> dv_x0(ncol);
    thrust::device_vector<float> dv_x(ncol);
    thrust::device_vector<float> dv_r0(ncol);
    thrust::device_vector<float> dv_w(ncol);

    thrust::device_vector<float> dv_H((m+1)*m);
    thrust::device_vector<float> dv_V((m+1)*ncol);
    thrust::device_vector<float> dv_y(m);
    
    thrust::device_ptr<float> dp_w = thrust::device_pointer_cast(dv_w.data());
    
    float *rp_x = thrust::raw_pointer_cast(dv_x.data());
    float *rp_x0 = thrust::raw_pointer_cast(dv_x0.data());
    float *rp_w = thrust::raw_pointer_cast(dv_w.data());
    float *rp_r0 = thrust::raw_pointer_cast(dv_r0.data());
    float *rp_H = thrust::raw_pointer_cast(dv_H.data());
    float *rp_V = thrust::raw_pointer_cast(dv_V.data());
    float *rp_y = thrust::raw_pointer_cast(dv_y.data());

    std::cout << "\nOur GMRES solution:\n\n";
    
    // get raw pointer
    csr_mat_t mat_A;
    mat_A.nrow = A.num_rows;
    mat_A.ncol = A.num_cols;
    mat_A.nnz = A.values.size();
    mat_A.value = thrust::raw_pointer_cast(A.values.data());
    mat_A.cindex = thrust::raw_pointer_cast(A.column_indices.data());
    mat_A.rowstart = thrust::raw_pointer_cast(A.row_offsets.data());

    vec_t vec_b;
    vec_b.value = thrust::raw_pointer_cast(b.data());
    vec_b.size = b.size();

    t_start = CycleTimer::currentSeconds();

    float b_norm2 = l2norm(b.begin(), b.end());

    thrust::fill(dv_x0.begin(), dv_x0.end(), .0);
    thrust::fill(dv_H.begin(), dv_H.end(), .0);

    while(nit < maxit) {

        gmres_compute_r0<<<blocks, threads>>>(rp_r0, mat_A, rp_x0, vec_b);
        float beta = l2norm(dv_r0.begin(), dv_r0.end());
        sax(1.0/beta, dv_r0.begin(), dv_r0.end(), dv_V.begin());

        innit = 0;

        for(int j = 0; j < m; j++){
            
            t_temp = CycleTimer::currentSeconds();

            s_mat_mul_x<<<blocks, threads>>>(rp_w, mat_A, (rp_V + j*ncol));

            for(size_t i = 0; i < j; i++){
                thrust::device_ptr<float> dp_vi = thrust::device_pointer_cast(rp_V+i*ncol);
                float dot = thrust::inner_product(dv_w.begin(), dv_w.end(), dp_vi, 0.f);
                float *hij = rp_H+i*m+j;

                sysax(dot, dp_vi, (dp_vi+ncol), dp_w);
                cudaMemcpy(hij, &dot, sizeof(float), cudaMemcpyHostToDevice);
            }
            float w_norm2 = l2norm(dv_w.begin(), dv_w.end());
            cudaMemcpy(rp_H+(j+1)*m+j, &w_norm2, sizeof(float), cudaMemcpyHostToDevice);


            thrust::device_ptr<float> dp_vjp1 = thrust::device_pointer_cast(rp_V+(j+1)*ncol);
            sax(1.0/w_norm2, dp_w, dp_w+ncol, dp_vjp1);

            t_krylov += CycleTimer::currentSeconds() - t_temp;
            t_temp = CycleTimer::currentSeconds();

            cudaMemcpy(H_host_data, rp_H, H_bytes, cudaMemcpyDeviceToHost);
            Matrix H_host_mat(m+1, m, H_host_data);
            Vector y_host_vec = leastSquareWithQR(H_host_mat, j+1, beta);
            cudaMemcpy(rp_y, y_host_vec.getData(), y_bytes, cudaMemcpyHostToDevice);
 
            t_lls += CycleTimer::currentSeconds() - t_temp;

            gmres_update_x<<<blocks, threads>>>(rp_x, rp_x0, rp_V, rp_y, j+1, ncol);
            gmres_compute_r0<<<blocks, threads>>>(rp_r0, mat_A, rp_x, vec_b);
            float res_norm = l2norm(dv_r0.begin(), dv_r0.end());

            nit++;
            innit++;
            if (res_norm < tol * b_norm2) {

                t_overall = CycleTimer::currentSeconds() - t_start;

                std::cout << "FGMRES converged to relative tolerance: "
                     << res_norm / b_norm2 << " at iteration " << nit
                     << " (out: " << outnit << ", in: " << innit << ")" << std::endl;

                sprintf(buf, "[%.3f] ms in Krylov \n", t_krylov * 1000);
                std::cout << buf;
                sprintf(buf, "[%.3f] ms in LLS \n", t_lls * 1000);
                std::cout << buf;
                sprintf(buf, "[%.3f] ms in Overall \n", t_overall * 1000);
                std::cout << buf;

                terminate = true;
                break;;
            }
        }
        if (terminate)
            break;

        thrust::copy(dv_x.begin(), dv_x.end(), dv_x0.begin());
        outnit++;
    }

    thrust::copy(dv_x.begin(), dv_x.end(), x.begin());

    // print out result for debug
    cusp::print(x);
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
    cusp::array1d<float, cusp::device_memory> x(A.num_cols, 0);     // 0
    cusp::array1d<float, cusp::device_memory> b(A.num_rows, 0);     // 1

    b[0] = 1;

    cusp::print(b);

    // get raw pointer
    // csr_mat_t csr_mat;
    // csr_mat.nrow = A.num_rows;
    // csr_mat.ncol = A.num_cols;
    // csr_mat.nnz = A.values.size();
    // csr_mat.value = thrust::raw_pointer_cast(A.values.data());
    // csr_mat.cindex = thrust::raw_pointer_cast(A.column_indices.data());
    // csr_mat.rowstart = thrust::raw_pointer_cast(A.row_offsets.data());

    // vec_t vec_x, vec_b;
    // vec_b.value = thrust::raw_pointer_cast(b.data());
    // vec_b.size = b.size();
    // vec_x.value = thrust::raw_pointer_cast(x.data());
    // vec_x.size = x.size();

    //gmres(csr_mat, vec_x, vec_b, 100, 1e-6, 1000);
    //cusp::print(x);

    gmres_new(A, x, b);
    gmres_ref(A, x, b);
}
