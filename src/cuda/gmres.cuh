#ifndef GMRES_CUDA_H__
#define GMRES_CUDA_H__

#include <cuda.h>

/* c struct for mat and vec */
typedef struct csr_mat_t {
    int *rowstart;
    int *cindex;
    float *value;
    int nrow;
    int ncol;
    int nnz;
} csr_mat_t;

typedef struct vec_t {
    float *value;
    int size;
} vec_t;

/* generic linear algebra kernel */
__global__ void s_mem_init(float *addr, float val, int N);
__global__ void s_x_sqrt(float *res, float *x, int N);
__global__ void s_x_div_a(float *res, float *x, float *a, int N);
__global__ void s_x_dot_y(float *res, float *x, float *y, int N);
__global__ void s_x_sub_ay(float *x, float *y, float *a, int N);
__global__ void s_mat_mul_x(float *res, csr_mat_t mat, float *x);


/* gmres related kernel */
__global__ void gmres_update_x(float *x, float *V, float *y, int m, int N);    
__global__ void gmres_compute_r0(float *r0, csr_mat_t mat, float *x, vec_t vec, float *beta);


/* host function prototype */
void cuda_handle_error(cudaError_t err, const char *file, int line);
void mem_log(float* device, int N);

#endif
