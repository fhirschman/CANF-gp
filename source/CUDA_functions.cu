#include "./../headers/data_structures.h"
#include "./../headers/CUDA_functions.cuh"

__host__
void freeAll_GPU(){
	cudaFree(POLY_COEFFS_GPU);
	cudaFree(POLY_EXPONS_GPU);
	cudaFree(iter_matrix_GPU);
	cudaFree(roots_matrix_GPU);
	cudaFree(cnumber_matrix_GPU);
	return;
}


__device__ double carg(const cuDoubleComplex& z) {
	return (double)atan2(cuCimag(z), cuCreal(z));
} 

__device__ cuDoubleComplex cpow(const cuDoubleComplex& z, const double& n) {
	return make_cuDoubleComplex((pow(cuCabs(z), n)*cos(n*carg(z))), (pow(cuCabs(z), n)*sin(n*carg(z))));
}

__device__ cuDoubleComplex clog_cu(const cuDoubleComplex& z) {
	return make_cuDoubleComplex( log(cuCabs(z)), carg(z) );
}

__device__ cuDoubleComplex cexp(const cuDoubleComplex& arg){
   cuDoubleComplex res;
   float s, c;
   float e = expf(arg.x);
   sincosf(arg.y, &s, &c);
   res.x = c * e;
   res.y = s * e;
   return res;
}


// Complex Gammafunction
__device__
cuDoubleComplex cgammln(double xx)
{
	cuDoubleComplex minusOne = make_cuDoubleComplex(-1.0f,0);
	cuDoubleComplex xxx,yyy,tmp,ser,res;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	xxx.x=xx;
	xxx.y=0;
	yyy.x=xx;
	yyy.y=0; 
	tmp.x = xxx.x + 5.5;
	tmp.y = 0;
	tmp=cuCadd(tmp,cuCmul(cuCmul(cuCadd(xxx,make_cuDoubleComplex(0.5,0)),make_cuDoubleComplex(log(tmp.x),0)),minusOne));
	ser=make_cuDoubleComplex(1.000000000190015,0);
	for (j=0;j<=5;j++){
		yyy = cuCadd(yyy,make_cuDoubleComplex(1.0f,0));
		ser =cuCadd(ser, cuCdiv( make_cuDoubleComplex(cof[j],0),yyy));
	}
	res.x = -tmp.x+cuCreal(clog_cu(make_cuDoubleComplex(2.5066282746310005*ser.x/xxx.x,0)));
	res.y = 0.0;
	return res;
}


__device__
cuDoubleComplex get_fdf_GPU(cuDoubleComplex& z, double* params, int n_terms, double* exponents,int FRAC_DER_MODE, double ALPHA_FRAC_DER ){
	cuDoubleComplex minusOne = make_cuDoubleComplex(-1.0f,0);
	cuDoubleComplex f = make_cuDoubleComplex(0,0);
	cuDoubleComplex df = make_cuDoubleComplex(0,0);
	cuDoubleComplex res = make_cuDoubleComplex(0,0);
	for(int i = n_terms-1; i>-1;i--){
		cuDoubleComplex rationalExponent = make_cuDoubleComplex(exponents[i],0);
		cuDoubleComplex temp1 = cpow(z,cuCreal(rationalExponent));
		f = cuCadd(f, cuCmul(make_cuDoubleComplex(params[i],0),temp1));		
	}
	// Normal derivative
	if(FRAC_DER_MODE==0){
		for(int i = n_terms-1; i>-1;i--){
			cuDoubleComplex rationalExponent = make_cuDoubleComplex(exponents[i],0);
			if(abs(params[i]) > 0.00000000001){
				cuDoubleComplex temp2 = cpow(z,cuCreal(cuCadd(rationalExponent,minusOne)));
				df = cuCadd(df,cuCmul(cuCmul(make_cuDoubleComplex(params[i],0),rationalExponent),temp2));
			}
		}
	}
	// Riemann-Liouville derivative
	else{
		for(int i = n_terms-1; i>-1;i--){
			if((abs(params[i]) > 0.00000000001)){
				if( (abs(exponents[i])  > 0.00000000001)) {
					cuDoubleComplex temp2 = cpow(z,abs(int(exponents[i]))-ALPHA_FRAC_DER);
					cuDoubleComplex temp3 = cgammln(abs(int(exponents[i]))-ALPHA_FRAC_DER+1.0);
					cuDoubleComplex temp4 = cgammln(abs(int(exponents[i]))+1.0);
					df = cuCadd(df, cuCmul(make_cuDoubleComplex( cuCreal( cuCdiv(cexp(temp4),cexp(temp3))),cuCimag(cuCdiv(cexp(temp4),cexp(temp3)))),temp2));
				}
				else{
					cuDoubleComplex temp2 = cpow(z,-ALPHA_FRAC_DER);
					cuDoubleComplex temp3 = cgammln(1.0-ALPHA_FRAC_DER);
					df = cuCadd(df,cuCdiv(temp2, make_cuDoubleComplex(cuCreal(cexp(temp3)),cuCimag(cexp(temp3)))));
				}
			}
		}
	}
	res = cuCmul(minusOne,cuCdiv(f,df));
	return res;
}

__global__
void solve_cnewton_iter_GPU(int* iter_matrix_GPU, cuDoubleComplex* cnumber_matrix_GPU, double accuracy, int max_iterations, double * params, int n_terms, double * rationalExponents,int N_X,int N_Y,int FRAC_DER_MODE, double ALPHA_FRAC_DER ){
	
	// printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
	int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
	int threadStride = blockDim.x * gridDim.x;
	cuDoubleComplex fdf = make_cuDoubleComplex(0,0);
	int matInd =threadIndex; 
	cuDoubleComplex z = make_cuDoubleComplex(cuCreal(cnumber_matrix_GPU[matInd]),cuCimag(cnumber_matrix_GPU[matInd]));
	for(int i = 0; i < max_iterations; i++){
		fdf = get_fdf_GPU(z,params,n_terms,rationalExponents,FRAC_DER_MODE,ALPHA_FRAC_DER);
		z = cuCadd(z,fdf);
		if(cuCabs(fdf) < accuracy){
			iter_matrix_GPU[matInd] = i;
			break;
		}
		else{
			iter_matrix_GPU[matInd] = max_iterations;
		}
	}
	return;
}

__global__
void solve_cnewton_roots_GPU(cuDoubleComplex* roots_matrix_GPU, cuDoubleComplex* cnumber_matrix_GPU, double accuracy, int max_iterations, double * params, int n_terms, double * rationalExponents,int N_X,int N_Y,int FRAC_DER_MODE, double ALPHA_FRAC_DER ){
	// printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
	int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
	int threadStride = blockDim.x * gridDim.x;
	cuDoubleComplex fdf = make_cuDoubleComplex(0,0);
	int matInd =threadIndex; 
	cuDoubleComplex z = make_cuDoubleComplex(cuCreal(cnumber_matrix_GPU[matInd]),cuCimag(cnumber_matrix_GPU[matInd]));
	for(int i = 0; i < max_iterations; i++){
		fdf = get_fdf_GPU(z,params,n_terms,rationalExponents,FRAC_DER_MODE,ALPHA_FRAC_DER );
		z = cuCadd(z,fdf);
		if(cuCabs(fdf) < accuracy){
			roots_matrix_GPU[matInd] = z;
			break;
		}
		else{
			roots_matrix_GPU[matInd] = make_cuDoubleComplex(0,0);
		}
	}
	return;
}


__host__
void setup_GPU(){
	numBlocks = (N_X*N_Y + BLOCK_SIZE - 1) / BLOCK_SIZE;
	POLY_COEFFS_GPU = new double[POLY_N_TERMS]();
	POLY_EXPONS_GPU = new double[POLY_N_TERMS]();
	cnumber_matrix_GPU = new cuDoubleComplex[N_X*N_Y];
	iter_matrix_GPU = new int[N_X*N_Y];
	roots_matrix_GPU = new cuDoubleComplex[N_X*N_Y];
	cudaMallocManaged(&POLY_COEFFS_GPU, (POLY_N_TERMS)*sizeof(double));
	cudaMallocManaged(&POLY_EXPONS_GPU, (POLY_N_TERMS)*sizeof(double));
	cudaMallocManaged(&iter_matrix_GPU, (N_X*N_Y)*sizeof(int));
	cudaMallocManaged(&cnumber_matrix_GPU, (N_X*N_Y)*sizeof(cuDoubleComplex));
	cudaMallocManaged(&roots_matrix_GPU, (N_X*N_Y)*sizeof(cuDoubleComplex));
	for (int i = 0; i < POLY_N_TERMS; i++){
		POLY_COEFFS_GPU[i] = POLY_COEFFS[i];
		POLY_EXPONS_GPU[i] = POLY_EXPONS[i];
	}
	return;
}

__host__
void fillCnumberMatrixGPU(){
	for (int i = 0; i<N_Y; i++){
		for (int j = 0; j<N_X; j++){
			cnumber_matrix_GPU[i*N_X+j] = make_cuDoubleComplex(real(cnumber_matrix[i][j]),(imag(cnumber_matrix[i][j])));
		}
	}
	return;
}

__host__
void getIterMatrixGPU(vector<vector<float>>& iterImage){
	solve_cnewton_iter_GPU<<<numBlocks,BLOCK_SIZE>>>(iter_matrix_GPU,cnumber_matrix_GPU,ACCURACY,N_ITER,POLY_COEFFS_GPU,POLY_N_TERMS,POLY_EXPONS_GPU,N_X,N_Y,FRAC_DER_MODE,ALPHA_FRAC_DER);
	cudaDeviceSynchronize();
	for(int i=0; i<N_Y; i++){
		vector<float> row;
		for(int j=0; j<N_X; j++){
			row.push_back(log(iter_matrix_GPU[j+i*N_X]));
		}
		iterImage.push_back(row);
	}
	return;
}

__host__
void getRootsMatrixGPU(vector<vector<float>>& rootsImage){
	solve_cnewton_roots_GPU<<<numBlocks,BLOCK_SIZE>>>(roots_matrix_GPU,cnumber_matrix_GPU,ACCURACY,N_ITER,POLY_COEFFS_GPU,POLY_N_TERMS,POLY_EXPONS_GPU,N_X,N_Y,FRAC_DER_MODE,ALPHA_FRAC_DER);
	cudaDeviceSynchronize();
	vector<complex<double>> roots;

	for(int i=0; i<N_Y; i++){
		vector<float> row;
		for(int j=0; j<N_X; j++){
			complex<double> root = {cuCreal(roots_matrix_GPU[j+i*N_X]),cuCimag(roots_matrix_GPU[j+i*N_X])};
			if(abs(real(root)) > 10*INTERVAL_LENGTH || abs(imag(root)) > 10*INTERVAL_LENGTH ){
				root = {0,0};
			} 
			if(i==0 && j==0){
				roots.push_back(root);
			}
			bool new_root = true;
			for (int k=0; k< roots.size(); k++){
				if( abs(real(root)-real(roots[k]))<0.01 && abs(imag(root)-imag(roots[k]))<0.01){
					roots_matrix[i][j]=k;
					row.push_back(roots_matrix[i][j]);
					new_root=false;
					break;
				}
			}
			if(new_root==true){
				roots_matrix[i][j]=roots.size()+1;
				row.push_back(roots_matrix[i][j]);
				roots.push_back(root);
			}
		}
		rootsImage.push_back(row);
	}
	cout << "Found " << roots.size() << " roots:" << endl;
	for (int k=0; k< roots.size(); k++){
		cout << " " << real(roots[k]) << " " << imag(roots[k]) << endl;
	}
	return;
}