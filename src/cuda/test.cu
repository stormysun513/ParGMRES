#include <iostream>

#include "gmres.cuh"

void debug(csr_mat_t mat){
        
    /* int dim = mat.ncol; */
 
    /* float *H; */
    /* float *V; */

    /* float *r0; */
    /* float *x; */
    /* float *w; */
    /* float *b; */
    /* float *tmp; */
    /* float *beta; */
    /* vec_t vec; */

    /* HANDLE_ERROR(cudaMalloc((void**)&H, (KRYLOV_M+1) * KRYLOV_M * sizeof(float))); */ 
    /* HANDLE_ERROR(cudaMalloc((void**)&V, (KRYLOV_M+1) * dim * sizeof(float))); */ 
    /* HANDLE_ERROR(cudaMalloc((void**)&r0, dim * sizeof(float))); */
    /* HANDLE_ERROR(cudaMalloc((void**)&x, dim * sizeof(float))); */
    /* HANDLE_ERROR(cudaMalloc((void**)&w, dim * sizeof(float))); */
    /* HANDLE_ERROR(cudaMalloc((void**)&b, dim * sizeof(float))); */
    /* HANDLE_ERROR(cudaMalloc((void**)&tmp, dim * sizeof(float))); */
    /* HANDLE_ERROR(cudaMalloc((void**)&beta, sizeof(float))); */

    /* int blocks = ROUND(dim, MAX_THREAD_NUM); */
    /* int threads = MAX_THREAD_NUM; */

    /* mem_init<<<blocks, threads>>>(x, 1.0, dim); */

    /* blocks = ROUND(dim, MAX_THREAD_NUM); */
    /* mem_init<<<blocks, threads>>>(b, .0, dim); */
 
    /* vec.value = b; */
    /* vec.size = dim; */

    /* HANDLE_ERROR(cudaDeviceSynchronize()); */

    /* compute_remainder<<<blocks, threads>>>(r0, mat, x, vec, tmp); */
    /* mem_log(r0, dim); */

    /* vector_sqrt<<<1,1>>>(beta, tmp, 1); */
    /* HANDLE_ERROR(cudaDeviceSynchronize()); */

    /* vector_divide_scalar<<<blocks, threads>>>(V, r0, beta, dim); */
    /* mem_log(V, dim); */
            

    /* int j = 0; */
    /* matrix_vector_multiply<<<blocks, threads>>>(w, mat, (V + j*dim)); */
    /* mem_log(w, dim); */
}

int main(){
    std::cout << "Test program\n";
    return 0;
}
