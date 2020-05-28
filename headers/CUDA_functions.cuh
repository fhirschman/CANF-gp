#pragma once
#include <cuComplex.h>

__host__ 
void freeAll_GPU();

__device__ 
cuDoubleComplex cgammln(double xx);

__device__ 
double carg(const cuDoubleComplex& z);

__device__ 
cuDoubleComplex clog_cu(const cuDoubleComplex& z);

__device__ 
cuDoubleComplex cpow(const cuDoubleComplex& z, const double& n);

__device__ 
cuDoubleComplex cexp (const cuDoubleComplex& arg);

__device__
cuDoubleComplex 
get_fdf_GPU(cuDoubleComplex& z, double* params, int n_terms, double* exponents,int FRAC_DER_MODE, double ALPHA_FRAC_DER);

__global__
void solve_cnewton_iter_GPU(int* iter_matrix_GPU, cuDoubleComplex* cnumber_matrix_GPU, double accuracy, int max_iterations, double* params, int n_terms, double * rationalExponents,int N_X,int N_Y,int FRAC_DER_MODE, double ALPHA_FRAC_DER);

__global__
void solve_cnewton_roots_GPU(cuDoubleComplex* roots_matrix_GPU, cuDoubleComplex* cnumber_matrix_GPU, double accuracy, int max_iterations, double* params, int n_terms, double* rationalExponents,int N_X,int N_Y,int FRAC_DER_MODE, double ALPHA_FRAC_DER);

__host__
void setup_GPU();

__host__
void fillCnumberMatrixGPU();

__host__
void getIterMatrixGPU(vector<vector<float>>& iterImage);

__host__
void getRootsMatrixGPU(vector<vector<float>>& rootsImage);